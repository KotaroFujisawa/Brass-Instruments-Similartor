"""Plot the digitized trumpet taper and the axial mesh used to resolve it."""
import os

from axisymmetric_mesh import refined_trumpet_mesh, trumpet_mesh
from trumpet import ROOT, body_profile, json, np

os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'tmp/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    out = ROOT/'output/trumpet_bell'
    os.makedirs(out, exist_ok=True)
    mouthpiece_end = json.loads((ROOT/'data/model.json').read_text())['cup_depth_m']+.061
    source_z, source_diameter = body_profile().T
    source_z = (source_z+mouthpiece_end)*1e3
    source_radius = source_diameter*500
    baseline = trumpet_mesh()
    refined = refined_trumpet_mesh()
    bell_z = refined['z'][-1]*1e3
    profiles = (
        (baseline, dict(color='#686d73', lw=1.3, ls='--', label='Baseline, linear interpolation')),
        (refined, dict(color='#00758a', lw=1.7, label='Refined, shape-preserving cubic')),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.0), layout='constrained')
    fig.suptitle('York B-flat trumpet bore and taper refinement', fontsize=16, fontweight='bold')

    lead = mouthpiece_end*1e3
    limits = ((0, bell_z+12), (8,325), (-120,2))
    shifts = (0,lead,bell_z)
    titles = ('Bore from mouthpiece to bell', 'Leadpipe taper', 'Final bell flare')
    xlabels = ('Axial distance from mouthpiece inlet (mm)',
               'Axial distance from mouthpipe inlet (mm)',
               'Position relative to bell rim (mm)')
    for ax, offset, xlim, title, label in zip(axes.flat,shifts,limits,titles,xlabels):
        for mesh, style in profiles:
            ax.plot(mesh['z']*1e3-offset,mesh['radius']*1e3,**style)
        select = (source_z-offset>=xlim[0]) & (source_z-offset<=xlim[1])
        ax.scatter(source_z[select]-offset,source_radius[select],s=15,facecolors='white',
                   edgecolors='#2b3940',linewidths=.8,zorder=3,label='Digitized source station')
        ax.set(xlim=xlim,title=title,xlabel=label,ylabel='Internal radius (mm)')
    axes[0,0].set_ylim(0,49)
    axes[0,0].legend(loc='upper left',frameon=False,fontsize=9)
    axes[0,1].set_ylim(4.6,5.6)
    axes[1,0].set_ylim(10,48)

    spacing = axes[1,1]
    for mesh, style in profiles:
        z = mesh['z']*1e3-bell_z
        spacing.stairs(np.diff(z),z,baseline=None,**style)
    spacing.set(xlim=(-120,2),ylim=(.9,3.8),title='Axial cell length near bell',
                xlabel='Position relative to bell rim (mm)',ylabel='Cell length (mm)')
    spacing.annotate('1.27 mm at the rim',xy=(-12,1.273),xytext=(-88,1.12),
                     arrowprops={'arrowstyle':'->','color':'#00758a'},color='#00758a')

    for ax in axes.flat:
        ax.grid(alpha=.22)
        ax.spines[['right','top']].set_visible(False)
    fig.savefig(out/'trumpet_taper_refinement.png',dpi=180)


if __name__=='__main__':
    main()
