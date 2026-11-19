// Axisymmetric Euler finite volumes. Density and total energy are stored
// relative to equilibrium, in units rho0 and rho0*c0^2, with velocity in c0.
struct Vec { float d, z, r, e; };
__device__ Vec operator+(Vec a, Vec b) { return {a.d+b.d,a.z+b.z,a.r+b.r,a.e+b.e}; }
__device__ Vec operator-(Vec a, Vec b) { return {a.d-b.d,a.z-b.z,a.r-b.r,a.e-b.e}; }
__device__ Vec operator*(float s, Vec a) { return {s*a.d,s*a.z,s*a.r,s*a.e}; }

struct Cell {
    int faces[6], neighbors[4];
    float sign_area_volume[6], radial_source, wall_loss, radius2;
};
struct Face { int left,right,kind; float nz,nr,eta_l,eta_r,area,radius2; };
struct Config {
    int nc,nf,inlet,bell,stride,mode,feedback,thermal_poles,radiation_order,reservoir_inflow;
    int radiation_modes,radiation_states;
    float dt,delta,area_in,area_bell,r_bell,guard,viscous_dc;
    float filters_b[5],filters_w[5];
    double lip_hz,force,rest_y,hold,attack,release;
    double lip_width,lip_thickness,lip_joint,lip_rest_x,lip_mass_product,lip_cup_area,lip_q;
    double attack_boost,attack_decay;
    double lip_gap_frequency_ratio,lip_gap_cubic_per_m2;
};
struct Clock { int step, error, fallback, pad; };
struct Boundary {
    double x,y,xpast,ypast,p,plip,ua,ul,old_ua,pm,in_h;
    float bell_memory,bell_predict,bell_derivative,bell_flux,bell_p;
    float bell_speed,bell_speed_predict,bell_speed_derivative;
};

constexpr float GAMMA=1.4f, P0=1.f/GAMMA, E0=P0/.4f;
constexpr double C0=343., RHO0=1.2, PUNIT=RHO0*C0*C0;

__device__ float pressure(Vec q) {
    return .4f*(q.e-.5f*(q.z*q.z+q.r*q.r)/(1.f+q.d));
}

__device__ float mc(float a,float b) {
    float s=fminf(.5f*fabsf(a+b),2.f*fminf(fabsf(a),fabsf(b)));
    return a*b>0 ? copysignf(s,a) : 0.f;
}

__device__ Vec slope(Vec a,Vec b) {
    return {mc(a.d,b.d),mc(a.z,b.z),mc(a.r,b.r),mc(a.e,b.e)};
}

__device__ Vec radial_slope(const Vec* values,const Cell& cell,Vec q) {
    int m=cell.neighbors[2],p=cell.neighbors[3];
    if (p<0) return {0,0,0,0};
    Vec prev=m<0 ? Vec{q.d,q.z,-q.r,q.e} : values[m];
    return slope(q-prev,values[p]-q);
}

__device__ Vec reconstruct(const Vec* values,const Cell* cells,int i,int dir,float side,float offset,float radius2) {
    const Cell& cell=cells[i];
    Vec q=values[i];
    if (dir==1) return q+side*radial_slope(values,cell,q);
    int m=cell.neighbors[0],p=cell.neighbors[1];
    Vec result=q;
    if (m>=0 && p>=0) {
        Vec s=slope(q-values[m],values[p]-q);
        float center=q.z*cell.radius2;
        s.z=mc(center-values[m].z*cells[m].radius2,values[p].z*cells[p].radius2-center)/cell.radius2;
        result=result+side*s;
    }
    if (offset!=0) result=result+offset*radial_slope(values,cell,q);
    // Reconstruct axial transport in the fitted rings, so a steady flux
    // through a steep taper does not appear as an O(rho*c*du) shock at a face.
    float z=result.z*cell.radius2/radius2;
    result.e+=.5f*(z*z-result.z*result.z)/(1.f+result.d);
    result.z=z;
    return result;
}

__device__ Vec physical_flux(Vec q,float p,float vn,float nz,float nr) {
    return {(1.f+q.d)*vn,q.z*vn+p*nz,q.r*vn+p*nr,(E0+P0+q.e+p)*vn};
}

__device__ Vec star_jump(Vec q,float p,float u,float side,float middle,float nz,float nr) {
    float rho=1.f+q.d, gap=middle-u;
    float drho=rho*gap/(side-middle);
    float dp=rho*(side-u)*gap;
    float de=(gap*(E0+P0+q.e+p)+middle*dp)/(side-middle);
    return {drho,drho*q.z/rho+(rho+drho)*gap*nz,
            drho*q.r/rho+(rho+drho)*gap*nr,de};
}

__device__ Vec riemann(Vec a,Vec b,Vec ca,Vec cb,float nz,float nr,Clock* clock) {
    float pa=pressure(a),pb=pressure(b);
    if (a.d<=-.99f || pa<=-.99f*P0) { a=ca; pa=pressure(a); atomicAdd(&clock->fallback,1); }
    if (b.d<=-.99f || pb<=-.99f*P0) { b=cb; pb=pressure(b); atomicAdd(&clock->fallback,1); }
    float ia=1.f/(1.f+a.d), ib=1.f/(1.f+b.d);
    float va=(nz*a.z+nr*a.r)*ia, vb=(nz*b.z+nr*b.r)*ib;
    float aa=sqrtf((1.f+GAMMA*pa)*ia),ab=sqrtf((1.f+GAMMA*pb)*ib);
    float sa=fminf(va-aa,vb-ab), sb=fmaxf(va+aa,vb+ab);
    if(sa>=0.f) return physical_flux(a,pa,va,nz,nr);
    if(sb<=0.f) return physical_flux(b,pb,vb,nz,nr);
    float la=(1.f+a.d)*(sa-va),lb=(1.f+b.d)*(sb-vb);
    float sm=(pb-pa+la*va-lb*vb)/(la-lb);
    float pstar=P0+pa+la*(sm-va);
    if(!(sm>sa && sm<sb && pstar>0.f)) {
        // HLLE fallback when the contact estimate leaves admissible states.
        return (1.f/(sb-sa))*(sb*physical_flux(a,pa,va,nz,nr)-sa*physical_flux(b,pb,vb,nz,nr)+sa*sb*(b-a));
    }
    return sm>=0.f ? physical_flux(a,pa,va,nz,nr)+sa*star_jump(a,pa,va,sa,sm,nz,nr)
                   : physical_flux(b,pb,vb,nz,nr)+sb*star_jump(b,pb,vb,sb,sm,nz,nr);
}

__device__ Vec wall(Vec q,float nz,float nr) {
    // Rusanov with the reflected state gives the exact zero mass/energy wall
    // flux, and a pressure traction with the wall-normal acoustic impedance.
    float inv=1.f/(1.f+q.d), p=pressure(q);
    float normal_momentum=nz*q.z+nr*q.r;
    float vn=normal_momentum*inv;
    float traction=p+normal_momentum*(vn+fabsf(vn)+sqrtf((1.f+GAMMA*p)*inv));
    return {0.f,nz*traction,nr*traction,0.f};
}

__device__ Vec open_face(Vec inside,float p,float u,bool inflow) {
    // The entering gas carries reference reservoir entropy, while exiting gas
    // carries local interior entropy. The pressure/flow characteristic is unchanged.
    float reference_p=inflow ? P0 : P0+pressure(inside);
    float reference_d=inflow ? 1.f : 1.f+inside.d;
    float density=reference_d*powf((P0+p)/reference_p,1.f/GAMMA);
    float v=inflow ? 0.f : inside.r/(1.f+inside.d);
    Vec face={density-1.f,density*u,density*v,p/.4f+.5f*density*(u*u+v*v)};
    return physical_flux(face,p,u,1.f,0.f);
}

__device__ void step_lips(Boundary* b,const Config* c,double t) {
    double w=c->lip_width, thickness=c->lip_thickness, joint=c->lip_joint, xeq=c->lip_rest_x;
    double envelope=1.;
    if(t<c->attack) envelope=.5-.5*cos(3.141592653589793*t/c->attack);
    else if(t>c->hold) envelope=.5+.5*cos(3.141592653589793*fmin(1.,(t-c->hold)/c->release));
    if(c->attack_boost>0. && t<c->attack+c->attack_decay) {
        double decay_phase=fmax(0.,t-c->attack)/c->attack_decay;
        envelope*=1.+.5*c->attack_boost*(1.+cos(3.141592653589793*decay_phase));
    }
    b->pm=c->force*envelope;
    double omega=6.283185307179586*c->lip_hz, dt=c->dt, inv_mass=omega*omega/(c->lip_mass_product*c->lip_hz);
    double stiff=omega*omega, damping=omega/(b->y>0. ? c->lip_q : .5);
    double gap_omega=omega*c->lip_gap_frequency_ratio;
    double gap_stiff=gap_omega*gap_omega, gap_damping=damping*c->lip_gap_frequency_ratio;
    double displacement=b->y-c->rest_y;
    double ax=2*w*inv_mass*(b->pm-b->p)*(joint-b->y)-stiff*(b->x-xeq);
    double ay=2*w*inv_mass*((b->pm-b->p)*b->x+thickness*b->plip)
             -gap_stiff*displacement*(1.+c->lip_gap_cubic_per_m2*displacement*displacement);
    if(b->y<0.) ay-=6*gap_stiff*b->y;
    double xn=(2*b->x-(1-damping*dt/2)*b->xpast+dt*dt*ax)/(1+damping*dt/2);
    double yn=(2*b->y-(1-gap_damping*dt/2)*b->ypast+dt*dt*ay)/(1+gap_damping*dt/2);
    b->ul=w*(b->x*(yn-b->ypast)-(b->y-joint)*(xn-b->xpast))/(2*dt);
    b->xpast=b->x; b->ypast=b->y; b->x=xn; b->y=yn;
    b->old_ua=b->ua;
}

__device__ float inlet_velocity(Boundary* b,const Config* c,float mean_h) {
    double w=c->lip_width, thickness=c->lip_thickness, scup=c->lip_cup_area;
    double h=c->feedback ? mean_h*PUNIT : 0.;
    double zb=RHO0*C0/c->area_in, area=2*w*b->y;
    double ua=0., gap=0.;
    if (area>1e-12) {
        double a=RHO0*(.5/(area*area)-1/(scup*area)+1/(scup*scup));
        double inertia=RHO0*thickness/(area*c->dt), l=inertia+zb;
        double drive=b->pm-h-zb*b->ul+inertia*b->old_ua;
        ua=2*drive/(l+sqrt(l*l+4*a*fabs(drive)));
        gap=-RHO0*ua*fabs(ua)*(1/(scup*area)-1/(scup*scup));
    }
    b->ua=ua; b->p=h+zb*(ua+b->ul); b->plip=b->p+gap; b->in_h=mean_h*PUNIT;
    return (ua+b->ul)/(c->area_in*C0);
}

__device__ float bell_reflection(Boundary* b,const Config* c,float outgoing,int stage) {
    float memory=stage ? b->bell_predict : b->bell_memory;
    if(c->radiation_order==2) {
        // Silva et al., causal unflanged reflection in scaled time c*t/a.
        float speed=stage ? b->bell_speed_predict : b->bell_speed;
        float acceleration=(outgoing-memory-1.393f*speed)/.457f;
        float step=c->delta/c->r_bell;
        if(stage) {
            b->bell_memory+=.5f*step*(speed+b->bell_derivative);
            b->bell_speed+=.5f*step*(acceleration+b->bell_speed_derivative);
        } else {
            b->bell_derivative=speed;
            b->bell_speed_derivative=acceleration;
            b->bell_predict=memory+step*speed;
            b->bell_speed_predict=speed+step*acceleration;
        }
        return -memory-.167f*speed;
    }
    float resistance=.61f*.61f/.25f, r_inf=(resistance-1)/(resistance+1);
    float rate=resistance/((resistance+1)*.61f*c->r_bell);
    float derivative=rate*(outgoing-memory);
    if(stage) b->bell_memory+=.5f*c->delta*(derivative+b->bell_derivative);
    else {
        b->bell_derivative=derivative;
        b->bell_predict=memory+c->delta*derivative;
    }
    return r_inf*outgoing-(1+r_inf)*memory;
}

__device__ float dot_product(const float* a,const float* b,int n) {
    float sum=0.f;
    for(int j=0;j<n;++j) sum=fmaf(a[j],b[j],sum);
    return sum;
}

__device__ void exterior_reflection(const Config* c,const float* coeff,float* state,
        const float* characteristic,float* incident,float* reflected,float* returned,int rings,int stage) {
    int lane=threadIdx.x, nm=c->radiation_modes, ns=c->radiation_states;
    const float* basis=coeff+nm*rings, *D=basis+nm*rings, *C=D+nm*nm;
    const float* Ad=C+nm*ns, *hold=Ad+ns*ns, *ramp=hold+ns*nm;
    float* next=state+ns, *first=next+ns;
    if(lane<nm) incident[lane]=.5f*dot_product(coeff+lane*rings,characteristic,rings);
    __syncthreads();
    if(lane<nm) {
        const float* at_stage=stage ? next : state;
        reflected[lane]=dot_product(D+lane*nm,incident,nm)+dot_product(C+lane*ns,at_stage,ns);
    }
    __syncthreads();
    if(stage) {
        if(lane<ns) {
            float correction=0.f;
            for(int j=0;j<nm;++j) correction=fmaf(ramp[lane*nm+j],incident[j]-first[j],correction);
            state[lane]=next[lane]+correction;
        }
    } else {
        if(lane<nm) first[lane]=incident[lane];
        if(lane<ns) {
            float predicted=0.f;
            for(int j=0;j<ns;++j) predicted=fmaf(Ad[j*ns+lane],state[j],predicted);
            next[lane]=predicted+dot_product(hold+lane*nm,incident,nm);
        }
    }
    if(lane<rings) returned[lane]=dot_product(basis+lane*nm,reflected,nm);
    __syncthreads();
}

extern "C" __global__ void fluxes(const Vec* current,const Vec* predicted,const Cell* cells,
        const Face* faces,Vec* flux,const Config* c,Clock* clock,Boundary* boundary,
        const float* radiation_coefficients,float* radiation_state,float* bell_trace,int record_bell,int stage) {
    const Vec* values=stage ? predicted : current;
    int internal_blocks=(c->inlet+blockDim.x-1)/blockDim.x;
    int lane=threadIdx.x;
    if (blockIdx.x<internal_blocks) {
        int i=blockIdx.x*blockDim.x+lane;
        if (i>=c->inlet) return;
        Face f=faces[i];
        Vec left=reconstruct(values,cells,f.left,f.kind>0,.5f,f.eta_l,f.radius2);
        if (f.kind==2) {
            if(left.d<=-.99f || pressure(left)<=-.99f*P0) left=values[f.left];
            flux[i]=wall(left,f.nz,f.nr);
        } else {
            Vec right=reconstruct(values,cells,f.right,f.kind,-.5f,f.eta_r,f.radius2);
            flux[i]=riemann(left,right,values[f.left],values[f.right],f.nz,f.nr,clock);
        }
        return;
    }
    bool inlet=blockIdx.x==internal_blocks;
    int start=inlet ? c->inlet : c->bell, stop=inlet ? c->bell : c->nf;
    __shared__ float radial_h[128], shared_condition;
    __shared__ float incident[16], reflected[16], returned[128];
    int i=start+lane;
    Face f;
    Vec q;
    if (i<stop) {
        f=faces[i]; q=values[inlet ? f.right : f.left];
        radial_h[lane]=pressure(q)+(inlet ? -1.f : 1.f)*q.z/(1.f+q.d);
    }
    __syncthreads();
    if(!inlet && c->radiation_modes>0 && (c->mode==0 || c->mode==3)) {
        int rings=stop-start;
        exterior_reflection(c,radiation_coefficients,radiation_state,radial_h,
                            incident,reflected,returned,rings,stage);
        // The mean bell observables are recorded by advance at the same
        // sampling steps. Reflection itself uses all rings every stage.
        if(lane==0 && stage && clock->step%c->stride==c->stride-1) {
            float flow=0.f, p=0.f;
            for(int j=0;j<rings;++j) {
                float a=.5f*radial_h[j], b=returned[j], area=faces[start+j].area;
                flow+=area*(a-b);
                p+=area*(a+b);
            }
            boundary->bell_flux=flow*C0;
            boundary->bell_p=p/c->area_bell*PUNIT;
        }
        if(i<stop) {
            float a=.5f*radial_h[lane], b=returned[lane], u=a-b;
            flux[i]=open_face(q,a+b,u,c->reservoir_inflow && u<0.f);
            if(record_bell && stage && clock->step%c->stride==c->stride-1)
                bell_trace[(clock->step/c->stride)*rings+lane]=u*C0;
        }
        return;
    }
    if(lane==0) {
        float integrated=0.f;
        for(int j=start;j<stop;++j) integrated+=faces[j].area*radial_h[j-start];
        if(inlet) {
            double t=(double)clock->step*c->dt;
            if(c->mode==0) {
                if(!stage) step_lips(boundary,c,t);
                shared_condition=inlet_velocity(boundary,c,integrated/c->area_in);
            } else if(c->mode==3) {
                double v=(t-.00045)/.00007;
                shared_condition=(c->force/PUNIT)*exp(-.5*v*v);
                boundary->in_h=integrated/c->area_in*PUNIT;
            }
        } else {
            // Reflection of the mean mode. Higher radial modes see a local
            // characteristic outlet; no external propagation grid is present.
            float outgoing=.5f*integrated/c->area_bell;
            shared_condition=bell_reflection(boundary,c,outgoing,stage);
            boundary->bell_flux=c->area_bell*(outgoing-shared_condition)*C0;
            boundary->bell_p=(outgoing+shared_condition)*PUNIT;
        }
    }
    __syncthreads();
    if(i>=stop) return;
    if(c->mode==1 || c->mode==2) flux[i]=inlet ? -1.f*wall(q,-1.f,0.f) : wall(q,1.f,0.f);
    else if(inlet) flux[i]=open_face(q,radial_h[lane]+shared_condition,shared_condition,
                                     c->reservoir_inflow && shared_condition>0.f);
    else {
        float u=.5f*radial_h[lane]-shared_condition;
        flux[i]=open_face(q,.5f*radial_h[lane]+shared_condition,u,c->reservoir_inflow && u<0.f);
        if(record_bell && stage && clock->step%c->stride==c->stride-1)
            bell_trace[(clock->step/c->stride)*(stop-start)+lane]=u*C0;
    }
}

__device__ Vec thermoviscous_loss(Vec q,float radius2,float4* memory,
                                  const float4* coefficients,const Config* c,int i) {
    float p=pressure(q), dc=c->viscous_dc/radius2;
    float z=dc*q.z, r=dc*q.r, thermal=0.f;
    for(int k=0;k<c->thermal_poles;++k) {
        int j=k*c->nc+i;
        float4 past=memory[j], coeff=coefficients[j];
        float dz=q.z-past.x, dr=q.r-past.y, dp=p-past.z;
        z+=coeff.z*dz; r+=coeff.z*dr; thermal+=coeff.w*dp;
        memory[j]=make_float4(past.x+coeff.x*dz,past.y+coeff.x*dr,past.z+coeff.y*dp,0.f);
    }
    q.z-=z; q.r-=r; q.e-=thermal/.4f;
    return q;
}

extern "C" __global__ void advance(Vec* current,Vec* predicted,const Cell* cells,const Vec* flux,
        void* loss_memory,const float4* loss_coefficients,const Config* c,Clock* clock,
        const Boundary* b,float* trace,int stage) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i>=c->nc) return;
    const Cell& cell=cells[i];
    Vec q=stage ? predicted[i] : current[i], rhs={0,0,0,0};
    for(int j=0;j<6;++j) {
        int id=cell.faces[j];
        if(id>=0) rhs=rhs+cell.sign_area_volume[j]*flux[id];
    }
    rhs.r-=cell.radial_source*pressure(q);
    Vec result=q-c->delta*rhs;
    if(stage) {
        result=.5f*(current[i]+result);
        if(c->thermal_poles>0 && cell.wall_loss>0.f) {
            result=thermoviscous_loss(result,cell.radius2,(float4*)loss_memory,loss_coefficients,c,i);
        } else if(cell.wall_loss>0.f) {
            float2* memory=(float2*)loss_memory;
            float z=0.f,r=0.f;
            for(int k=0;k<5;++k) {
                float2 past=memory[k*c->nc+i];
                float dz=result.z-past.x, dr=result.r-past.y;
                z+=c->filters_w[k]*dz; r+=c->filters_w[k]*dr;
                memory[k*c->nc+i]=make_float2(past.x+c->filters_b[k]*dz,past.y+c->filters_b[k]*dr);
            }
            result.z-=cell.wall_loss*z; result.r-=cell.wall_loss*r;
        }
        current[i]=result;
        float inv=1.f/(1.f+result.d), sound2=(1.f+GAMMA*pressure(result))*inv;
        if(!(inv>0.f && sound2>0.f && isfinite(sound2))) atomicOr(&clock->error,1);
        if((fabsf(result.z)+fabsf(result.r))*inv+sqrtf(sound2)>c->guard) atomicOr(&clock->error,2);
        if(i==0) {
            if(clock->step%c->stride==c->stride-1) {
                float* row=trace+11*(clock->step/c->stride);
                row[0]=b->p; row[1]=b->ua+b->ul; row[2]=b->x; row[3]=b->y;
                row[4]=b->ua; row[5]=b->ul; row[6]=b->pm; row[7]=b->plip;
                row[8]=b->in_h; row[9]=b->bell_flux; row[10]=b->bell_p;
            }
            ++clock->step;
        }
    } else predicted[i]=result;
}
