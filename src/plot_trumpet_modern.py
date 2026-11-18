"""Plot input impedance after running verify_trumpet_modern.py."""
import os

from render_trumpet_modern import OUT, ROOT, np

os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'tmp/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    with np.load(OUT/'input_impedance.npz') as data:
        frequency = data['frequency_hz']
        previous = abs(data['previous_pa_s_m3'])/1e6
        updated = abs(data['modern_pa_s_m3'])/1e6

    fig, ax = plt.subplots(figsize=(8.4, 3.4), layout='constrained')
    ax.plot(frequency, previous, ls='--', lw=1.25, color='#777777',
            label='Previous wall and bell models')
    ax.plot(frequency, updated, lw=1.45, color='#1874a4',
            label='Updated thermoviscous wall and bell models')
    ax.set(xlim=(50, 1500), ylim=(0, 120), xlabel='Frequency (Hz)',
           ylabel='Input impedance magnitude (MPa·s/m³)',
           title='Trumpet input impedance, open fingering')
    ax.legend(frameon=False)
    ax.grid(alpha=.22)
    fig.savefig(OUT/'trumpet_modern_impedance.png', dpi=180)
    plt.close(fig)


if __name__=='__main__':
    main()
