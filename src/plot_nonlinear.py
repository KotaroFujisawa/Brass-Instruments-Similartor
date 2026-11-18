"""Comparison figure for the fixed-source Burgers extension."""
from nonlinear_propagation import *
import os
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    folder=ROOT/'output/nonlinear'
    report=json.loads((folder/'measurements.json').read_text())
    fields=np.load(folder/'Bb4_fields.npz')
    n=fields['nonlinear_radiated_pa'];l=fields['linear_radiated_pa'];fs=48000
    time=np.arange(len(n))/fs
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
    fig.suptitle('TRUMPET / NONLINEAR BORE PROPAGATION\nGeneralised Burgers equation - prescribed incident wave, matched linear reference',fontsize=16)
    select=(time>=.34)&(time<=.347)
    axs[0,0].plot((time[select]-.34)*1000,l[select],label='Linear',color='#6c8092',lw=1.7)
    axs[0,0].plot((time[select]-.34)*1000,n[select],label='Nonlinear',color='#c56f1d',lw=1.7)
    axs[0,0].set(xlabel='Time (ms)',ylabel='Radiated pressure (Pa)',title='Bb4: changed waveform and steeper pressure fronts')
    axs[0,0].legend()
    seg=slice(int(.3*fs),int(.56*fs))
    fq=rfftfreq(131072,1/fs)
    spectra=[]
    for x in [l,n]:spectra.append(abs(rfft(x[seg]*np.hanning(len(x[seg])),131072)))
    ref=spectra[0].max()
    for s,label,col in zip(spectra,['Linear','Nonlinear'],['#6c8092','#c56f1d']):
        axs[0,1].plot(fq,20*np.log10(np.maximum(s/ref,1e-8)),label=label,color=col,lw=1.)
    axs[0,1].set(xlim=(0,10000),ylim=(-90,12),xlabel='Frequency (Hz)',ylabel='Level vs linear peak (dB)',title='Bb4: harmonic enrichment (same physical gain)')
    axs[0,1].legend()
    notes=report['notes'];idx=np.arange(len(notes))
    axs[1,0].bar(idx,[x['hf_fraction_ratio_db'] for x in notes],color='#c56f1d',width=.65)
    axs[1,0].axhline(0,color='gray',lw=.8)
    axs[1,0].set(xticks=idx,xticklabels=[x['note'] for x in notes],ylabel='High-frequency fraction ratio (dB)',title='Increase in energy fraction above 3 kHz')
    axs[1,0].tick_params(axis='x',rotation=50)
    axs[1,1].plot(idx,[x['cents'] for x in notes],'o-',color='#315e83')
    axs[1,1].axhline(0,color='gray',lw=.8)
    axs[1,1].axhline(5,color='gray',lw=.7,ls=':');axs[1,1].axhline(-5,color='gray',lw=.7,ls=':')
    axs[1,1].set(xticks=idx,xticklabels=[x['note'] for x in notes],ylabel='Pitch error (cents)',title='Measured radiated pitch / A4 = 440 Hz')
    axs[1,1].tick_params(axis='x',rotation=50)
    fig.savefig(folder/'nonlinear_comparison.png',dpi=160)
    plt.close(fig)
    print('Saved nonlinear_comparison.png')


if __name__=='__main__':main()
