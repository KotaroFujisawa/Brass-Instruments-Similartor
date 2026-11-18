"""Pressure field and validation plots for the axisymmetric model."""
from render_axisymmetric import OUT
from verify_axisymmetric import pressure
from trumpet import C, RHO, ROOT, json, np
import os

os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


def main():
    field=np.load(OUT/'Bb4_field.npz');q=field['field']
    z,radii,rings=field['z'],field['radius'],field['rings']
    p=RHO*C*C*pressure(q)
    tests=np.load(OUT/'solver_curves.npz')
    stats=json.loads((OUT/'measurements.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(12,7.6),constrained_layout=True)
    fig.suptitle('Axisymmetric Euler trumpet / concert Bb3 to Bb5',fontsize=15)
    ax=axes[0,0];polys=[];values=[]
    ids=field['offsets']
    for i,n in enumerate(rings):
        if z[i+1]<z[-1]-.19:continue
        for j in range(n):
            a,b=j/n,(j+1)/n
            polys.append([(z[i],a*radii[i]),(z[i+1],a*radii[i+1]),
                          (z[i+1],b*radii[i+1]),(z[i],b*radii[i])])
            values.append(p[ids[i]+j])
    collection=PolyCollection(polys,array=np.array(values),cmap='RdBu_r',linewidths=0.)
    bound=max(abs(np.array(values)));collection.set_clim(-bound,bound)
    ax.add_collection(collection);ax.autoscale_view();fig.colorbar(collection,ax=ax,label='Gauge pressure (Pa)')
    ax.set(xlabel='Axial position (m)',ylabel='Radius (m)',title='Bb4: pressure in the bell at 0.41 s')
    ax=axes[0,1];rr=np.diff(ids)
    ax.step(.5*(z[1:]+z[:-1]),rr,where='mid',color='#355b83')
    ax.set(xlabel='Axial position (m)',ylabel='Annular cells per section',title='Radial mesh through the mouthpiece and bore')
    ax=axes[1,0]
    ax.plot(tests['x_exact'],tests['p_exact'],color='#282828',label='Sod exact',lw=1.3)
    ax.plot(tests['x_coarse'],tests['p_coarse'],color='#8799a6',label='200 axial cells')
    ax.plot(tests['x_fine'],tests['p_fine'],color='#c26521',ls='--',label='400 axial cells')
    ax.set(xlabel='Position (m)',ylabel='Dimensionless pressure',title='Shock tube: rarefaction and shock');ax.legend(fontsize=9)
    ax=axes[1,1];notes=stats['notes'];ids=np.arange(len(notes))
    ax.axhline(0,color='#6d7682',lw=.7);ax.plot(ids,[s['cents'] for s in notes],'o-',color='#355b83')
    ax.set(xticks=ids,xticklabels=[s['note'] for s in notes],ylabel='Pitch error (cents)',title='Forte scale: radiated pitch');ax.tick_params(axis='x',rotation=45)
    fig.savefig(OUT/'axisymmetric_diagnostics.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
