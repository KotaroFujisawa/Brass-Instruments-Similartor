"""Plot the reproduced bore, impedance, pressure snapshot, and musical pitch."""
import os

from trombone_axisymmetric import OUT, profile
from trumpet import C, RHO, ROOT, json, np
from verify_axisymmetric import pressure

os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


def main():
    data=json.loads((OUT/'measurements.json').read_text())
    with np.load(OUT/'impedance_curves.npz') as a:
        f,z1,z7=a['f_hz'],a['pos1'],a['pos7']
    with np.load(OUT/'Bb3_field.npz') as a:
        z,r,counts,offsets=a['z'],a['radius'],a['rings'],a['offsets']
        p=pressure(a['state'])*RHO*C*C
    fig,axs=plt.subplots(2,2,figsize=(12,7.6),constrained_layout=True)
    fig.suptitle('Conn 8H tenor trombone + Denis Wick 5AL / axisymmetric Euler',fontsize=14)
    ax=axs[0,0]
    for position,color in [(1,'#215d84'),(7,'#bd5b28')]:
        x,radius=profile(position)
        ax.plot(x,1000*radius,color=color,label=f'Position {position}')
    ax.set(xlabel='Distance along unfolded bore (m)',ylabel='Bore radius (mm)',title='Reconstructed bore and slide')
    ax.legend();ax.set_ylim(0,116)
    ax=axs[0,1]
    ax.plot(f,z1/1e6,color='#215d84',label='Position 1')
    ax.plot(f,z7/1e6,color='#bd5b28',label='Position 7',alpha=.85)
    ax.set(xlim=(10,650),xlabel='Frequency (Hz)',ylabel='Input impedance (MPa s/m³)',title='Small-signal impedance')
    ax.legend()
    ax=axs[1,0];boxes=[];p_values=[]
    for k,n in enumerate(counts):
        if z[k+1]<z[-1]-.32:continue
        for j in range(n):
            a,b=j/n,(j+1)/n
            boxes.append([(z[k],a*r[k]),(z[k+1],a*r[k+1]),(z[k+1],b*r[k+1]),(z[k],b*r[k])])
            p_values.append(p[offsets[k]+j])
    colored=PolyCollection(boxes,array=np.array(p_values),cmap='RdBu_r',linewidths=0.)
    limit=max(abs(np.array(p_values)));colored.set_clim(-limit,limit)
    ax.add_collection(colored);ax.autoscale_view();fig.colorbar(colored,ax=ax,label='Gauge pressure (Pa)')
    ax.set(xlabel='Axial distance (m)',ylabel='Radius (m)',title='Bb3 pressure in the bell at 0.42 s')
    ax=axs[1,1];notes=data['notes'];j=np.arange(len(notes))
    ax.axhline(0,color='#677180',lw=.8)
    ax.plot(j,[note['cents'] for note in notes],'o-',color='#215d84')
    ax.set(xticks=j,xticklabels=[note['note'] for note in notes],ylabel='Pitch error (cents)',title='Measured pitch of the radiated scale')
    ax.tick_params(axis='x',rotation=45)
    fig.savefig(OUT/'trombone_diagnostics.png',dpi=160)
    plt.close(fig)


if __name__=='__main__':main()
