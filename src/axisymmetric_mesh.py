"""Annular finite volumes with conical radial faces and nested rings."""
from trumpet import C, ROOT, geometry, json, np, optimize

CELL = np.dtype([('faces', 'i4', (6,)), ('neighbors', 'i4', (4,)),
                 ('sign_area_volume', 'f4', (6,)), ('radial_source', 'f4'),
                 ('wall_loss', 'f4'), ('radius2', 'f4')], align=True)
FACE = np.dtype([('left', 'i4'), ('right', 'i4'), ('kind', 'i4'),
                 ('nz', 'f4'), ('nr', 'f4'), ('eta_l', 'f4'), ('eta_r', 'f4'),
                 ('area', 'f4'), ('radius2', 'f4')], align=True)


def trumpet_mesh(valves=(), dz=.0035, dr=.003, tuning=0.):
    depth = json.loads((ROOT/'data/model.json').read_text())['cup_depth_m']
    _, _, x, radius = geometry(valves, depth, tuning)
    # Include the throat endpoints and retain finer axial resolution in the cup.
    anchors = [0., depth, depth+.006, depth+.061, x[-1]-.11, x[-1]]
    spacing = [min(.0018,dz), min(.002,dz), dz, dz, min(.0025,dz)]
    z = np.concatenate([np.linspace(a,b,1+int(np.ceil((b-a)/h)))[:-1]
                        for a,b,h in zip(anchors[:-1],anchors[1:],spacing)] + [x[-1:]])
    return annular_mesh(z, np.interp(z,x,radius), dr)


def annular_mesh(z, radius, dr=.003, loss=True, min_rings=2):
    z, radius = np.asarray(z), np.asarray(radius)
    widths = np.diff(z)
    center_r = (radius[:-1]+radius[1:])/2
    rings = (2**np.maximum(int(np.log2(min_rings)),np.round(np.log2(center_r/dr)))).astype(int)
    # Adjacent columns differ by at most one ring bisection.
    for i in range(1,len(rings)):
        rings[i] = max(rings[i], rings[i-1]//2)
    for i in range(len(rings)-2,-1,-1):
        rings[i] = max(rings[i], rings[i+1]//2)
    offsets = np.r_[0,np.cumsum(rings)]
    cells = np.zeros(offsets[-1],dtype=CELL)
    cells['faces'] = -1
    cells['neighbors'] = -1
    volumes = np.empty(len(cells)); centers = np.empty((len(cells),2))
    face_count = np.zeros(len(cells),dtype=int)
    radial_source = np.zeros(len(cells)); total_area = np.zeros(len(cells))
    faces = []

    def add(left,right,kind,nz,nr,area,eta_l=0.,eta_r=0.,radius2=0.):
        fid = len(faces)
        faces.append((left,right,kind,nz,nr,eta_l,eta_r,area,radius2))
        for index,sign in [(left,1),(right,-1)]:
            if index < 0:
                continue
            slot = face_count[index]; face_count[index] += 1
            cells['faces'][index,slot] = fid
            cells['sign_area_volume'][index,slot] = sign*area/volumes[index]
            radial_source[index] += sign*area*nr/volumes[index]
            total_area[index] += area

    for iz,n in enumerate(rings):
        a,b = radius[iz:iz+2]
        full_volume = np.pi*widths[iz]*(a*a+a*b+b*b)/3
        for j in range(n):
            cid = offsets[iz]+j
            eta = (j+.5)/n
            volumes[cid] = full_volume*(2*j+1)/n**2
            centers[cid] = [(z[iz]+z[iz+1])/2, eta*center_r[iz]]
            left = offsets[iz-1]+int(eta*rings[iz-1]) if iz>0 else -1
            right = offsets[iz+1]+int(eta*rings[iz+1]) if iz+1<len(rings) else -1
            cells['neighbors'][cid] = [left,right,cid-1 if j else -1,cid+1 if j+1<n else -1]
            cells['radius2'][cid] = center_r[iz]**2
            if loss:
                cells['wall_loss'][cid] = 4*np.sqrt(np.pi)*1.045/(632.8*center_r[iz])

    for iz,n in enumerate(rings):
        a,b = radius[iz:iz+2]
        for j in range(1,n+1):
            eta = j/n
            slope = eta*(b-a)/widths[iz]
            norm = np.hypot(1,slope)
            area = 2*np.pi*eta*center_r[iz]*widths[iz]*norm
            add(offsets[iz]+j-1,offsets[iz]+j if j<n else -1,
                1 if j<n else 2,-slope/norm,1/norm,area)
        if iz+1 == len(rings):
            continue
        m = rings[iz+1]; sections = max(m,n)
        for j in range(sections):
            eta = (j+.5)/sections
            l,r = int(eta*n),int(eta*m)
            area = np.pi*b*b*(2*j+1)/sections**2
            add(offsets[iz]+l,offsets[iz+1]+r,0,1.,0.,area,eta*n-l-.5,eta*m-r-.5,b*b)
    inlet = len(faces)
    for j in range(rings[0]):
        n = rings[0]
        add(-1,j,3,1.,0.,np.pi*radius[0]**2*(2*j+1)/n**2)
    bell = len(faces)
    for j in range(rings[-1]):
        n = rings[-1]
        add(offsets[-2]+j,-1,4,1.,0.,np.pi*radius[-1]**2*(2*j+1)/n**2)
    cells['radial_source'] = radial_source
    # Deflate the eigenvalue-based finite-volume limit to leave room for
    # finite-amplitude wave speed, MUSCL, the axis, and nonorthogonal faces.
    length_scale = 2*np.min(volumes/total_area)
    return dict(cells=cells,faces=np.array(faces,dtype=FACE),inlet=inlet,bell=bell,
                z=z,radius=radius,rings=rings,offsets=offsets,centers=centers,
                volumes=volumes,length_scale=float(length_scale))


def loss_fit():
    theta = 2*np.pi*np.geomspace(18,90000,5)
    s = 2j*np.pi*np.geomspace(125,18000,180)
    a = s[:,None]/((s[:,None]+theta)*np.sqrt(s[:,None]))
    w,_ = optimize.nnls(np.r_[a.real,a.imag],np.r_[np.ones(len(s)),np.zeros(len(s))])
    return theta,w,float(np.max(abs(a@w-1)))
