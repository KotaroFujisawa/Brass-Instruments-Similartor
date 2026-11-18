"""Scientific diagnostics for two-way nonlinear feedback."""
from bidirectional import *
import os
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    folder=ROOT/'output/bidirectional'
    report=json.loads((folder/'measurements.json').read_text())
    fields=np.load(folder/'Bb4_coupled.npz');fs=int(fields['fs'])
    n=fields['nonlinear'];l=fields['linear']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(3,2,figsize=(13,12),constrained_layout=True)
    fig.suptitle('TRUMPET / TWO-WAY NONLINEAR FEEDBACK\nSelf-excited 2D lips + mouthpiece + bore + causal reflected waves',fontsize=16)
    colors=['#607c91','#c36620'];labels=['Linear bore (beta=0)','Nonlinear bore (beta=1.2)']
    # Use identical absolute times; resulting phase shift is part of feedback.
    section=slice(int(.34*fs),int(.347*fs));t=np.arange(len(n))/fs
    for tr,color,label in zip([l,n],colors,labels):
        axs[0,0].plot((t[section]-.34)*1000,tr[section,8]/1000,color=color,label=label,lw=1.4)
        axs[0,1].plot(tr[section,2]*1000,tr[section,3]*1000,color=color,label=label,lw=1.4)
    axs[0,0].set(xlabel='Time (ms)',ylabel='Incoming pressure p- (kPa)',title='Bb4: actual nonlinear return reaching the lips')
    axs[0,0].legend(fontsize=9)
    axs[0,1].set(xlabel='Lip x (mm)',ylabel='Lip y (mm)',title='Bb4: feedback changes the lip orbit')
    axs[0,1].axhline(0,color='gray',ls=':',lw=.8)
    afs=48000;ts=np.arange(len(fields['nonlinear_radiated_pa']))/afs
    sel=(ts>=.34)&(ts<=.347);seg=slice(int(.30*afs),int(.56*afs));fq=np.fft.rfftfreq(131072,1/afs)
    spectra=[]
    for key,col,label in zip(['linear_radiated_pa','nonlinear_radiated_pa'],colors,labels):
        y=fields[key]
        axs[1,0].plot((ts[sel]-.34)*1000,y[sel],color=col,label=label,lw=1.4)
        spectra.append(abs(np.fft.rfft(y[seg]*np.hanning(len(y[seg])),131072)))
    reference=max(spectra[0])
    for spec,col,label in zip(spectra,colors,labels):
        axs[1,1].plot(fq,20*np.log10(np.maximum(spec/reference,1e-8)),color=col,label=label,lw=1)
    axs[1,0].set(xlabel='Time (ms)',ylabel='Radiated pressure (Pa)',title='Bb4: bell radiation at 1 m (monopole approximation)')
    axs[1,1].set(xlim=(0,10000),ylim=(-90,10),xlabel='Frequency (Hz)',ylabel='Level vs linear peak (dB)',title='Bb4: radiated spectrum, identical physical gain')
    notes=report['notes'];idx=np.arange(len(notes))
    ratios=[10*np.log10(x['nonlinear']['hf_fraction_above_3khz']/x['linear']['hf_fraction_above_3khz']) for x in notes]
    axs[2,0].bar(idx,ratios,color=colors[1],width=.65)
    axs[2,0].axhline(0,color='gray',lw=.8)
    axs[2,0].set(xticks=idx,xticklabels=[x['note'] for x in notes],ylabel='HF fraction ratio (dB)',title='Change in energy fraction above 3 kHz')
    axs[2,1].plot(idx,[x['cents'] for x in notes],'o-',color='#315e83')
    axs[2,1].axhline(0,color='gray',lw=.8)
    axs[2,1].axhline(5,color='gray',ls=':',lw=.8);axs[2,1].axhline(-5,color='gray',ls=':',lw=.8)
    axs[2,1].set(xticks=idx,xticklabels=[x['note'] for x in notes],ylabel='Pitch error (cents)',title='Delivered nonlinear scale / A4 = 440 Hz')
    for ax in axs[2,:]:ax.tick_params(axis='x',rotation=50)
    fig.savefig(folder/'diagnostics.png',dpi=160);plt.close(fig)
    print('Saved diagnostics.png',flush=True)


if __name__=='__main__':main()
