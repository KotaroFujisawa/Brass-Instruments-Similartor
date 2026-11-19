"""Measure, plot, and export the trumpet bell comparison from pressure traces."""
import json
import os
from pathlib import Path

from render_axisymmetric import ROOT, TRACE_FS
from render_bidirectional import measure_audio
from trumpet import measure, np
from scipy.io import wavfile

AUDIO_FS = 48000
NOTES = ('Bb4', 'Bb5', 'Bb6')
VARIANTS = ('previous', 'fine_taper', 'exterior')
LABELS = ('Original taper and bell', 'Refined taper, original bell', 'Refined taper, exterior bell')
COLORS = ('#53667a', '#a45b14', '#047e7b')
P_REF = 20e-6


def rms(values):
    return float(np.std(values, dtype=np.float64))


def harmonic_rms(wave, fs, pitch, limit=10):
    """Estimate individual line amplitudes with a Hann-windowed spectrum."""
    if not pitch:
        return []
    n = len(wave)
    win = np.hanning(n)
    nfft = 1 << (4*n-1).bit_length()
    f = np.fft.rfftfreq(nfft, 1/fs)
    rms_spectrum = np.sqrt(2)*abs(np.fft.rfft((wave-wave.mean())*win, nfft))/win.sum()
    result = []
    for order in range(1, limit+1):
        center = order*pitch
        if center >= min(20000, .49*fs):
            break
        half = max(2*fs/n, .005*center)
        lo, hi = np.searchsorted(f, [center-half, center+half])
        if hi <= lo:
            continue
        index = lo+np.argmax(rms_spectrum[lo:hi])
        amplitude = float(rms_spectrum[index])
        frequency = float(f[index])
        neighbors = rms_spectrum[index-1:index+2]
        if len(neighbors)==3 and np.all(neighbors>0) and amplitude>max(neighbors[0],neighbors[2]):
            left, center, right = np.log(neighbors)
            delta = .5*(left-right)/(left-2*center+right)
            amplitude = float(np.exp(center-.25*(left-right)*delta))
            frequency += float(delta*fs/nfft)
        result.append(dict(order=order, frequency_hz=frequency, rms_pa=amplitude,
                           spl_db=float(20*np.log10(amplitude/P_REF)) if amplitude>0 else None))
    return result


def analyze_record(note, variant, record):
    trace = np.asarray(record['trace'])
    audio = {int(angle): np.asarray(wave) for angle, wave in record['audio'].items()}
    if variant!='exterior':
        # A monopole is angle-independent, so one waveform serves both views.
        audio.setdefault(45, audio[0])
    duration = float(record.get('duration_s', len(audio[0])/AUDIO_FS))
    start, end = map(float, record.get('window_s', (max(.04, duration-.36), duration-.16)))
    if (trace.ndim!=2 or trace.shape[1]<10 or not np.isfinite(trace).all()
            or len(trace)!=round(duration*TRACE_FS) or not 0<=start<end<=duration):
        raise ValueError(f'Invalid trace or analysis window for {note} / {variant}.')
    if (not {0,45}.issubset(audio) or
            any(w.ndim!=1 or len(w)!=round(duration*AUDIO_FS) or not np.isfinite(w).all()
                for w in audio.values())):
        raise ValueError(f'Invalid 0°/45° audio for {note} / {variant}.')
    target = float(record['target_hz'])
    cup = trace[round(start*TRACE_FS):round(end*TRACE_FS)]
    p_cup = measure(cup, fs=TRACE_FS, start=0., max_hz=2450)
    p_audio = measure_audio(audio[0],start,end,max_hz=2450)
    # Use one fundamental for all angles. This keeps weak off-axis harmonics
    # registered to the same acoustic source even if a pitch estimator skips them.
    pitch = p_cup['frequency'] or p_audio['frequency']
    radiation = {}
    for angle, wave in sorted(audio.items()):
        steady = wave[round(start*AUDIO_FS):round(end*AUDIO_FS)]
        level = rms(steady)
        radiation[str(angle)] = dict(rms_pa=level, spl_db=float(20*np.log10(level/P_REF)) if level>0 else None,
                                    harmonics=harmonic_rms(steady,AUDIO_FS,pitch))
    h0 = radiation['0']['harmonics']
    h45 = radiation['45']['harmonics']
    for h, normal in zip(h45,h0):
        h['relative_to_0deg_db'] = (h['spl_db']-normal['spl_db']
                                   if h['spl_db'] is not None and normal['spl_db'] is not None else None)
    row = dict(
        note=note, variant=variant, target_hz=target, duration_s=duration, window_s=[start,end],
        cup_frequency_hz=float(p_cup['frequency']) or None, cup_periodicity=float(p_cup['periodicity']),
        cup_rms_pa=rms(cup[:,0]), cup_peak_to_peak_pa=float(np.ptp(cup[:,0])),
        bell_volumetric_flow_rms_m3_s=rms(cup[:,9]),
        frequency_hz=float(p_audio['frequency']) or None, periodicity=float(p_audio['periodicity']),
        pitch_cents=float(1200*np.log2(p_audio['frequency']/target)) if p_audio['frequency'] else None,
        radiation=radiation, diagnostics=record.get('diagnostics', {}),
    )
    return row, audio


def pitch_delta(first, second, field):
    if first[field] and second[field]:
        return float(1200*np.log2(second[field]/first[field]))
    return None


def comparisons(rows):
    result = {}
    for note in NOTES:
        previous, fine, exterior = (rows[note,variant] for variant in VARIANTS)
        result[note] = dict(
            fine_taper_minus_previous_pitch_cents=pitch_delta(previous,fine,'frequency_hz'),
            exterior_minus_fine_taper_pitch_cents=pitch_delta(fine,exterior,'frequency_hz'),
            fine_taper_minus_previous_cup_pitch_cents=pitch_delta(previous,fine,'cup_frequency_hz'),
            exterior_minus_fine_taper_cup_pitch_cents=pitch_delta(fine,exterior,'cup_frequency_hz'),
            fine_taper_minus_previous_cup_rms_pa=fine['cup_rms_pa']-previous['cup_rms_pa'],
            exterior_minus_fine_taper_cup_rms_pa=exterior['cup_rms_pa']-fine['cup_rms_pa'],
        )
    return result


def plot_results(rows, out):
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/mplconfig'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1,3,figsize=(12,4),layout='constrained',sharey=True)
    fig.suptitle('Radiated harmonic levels at 1 m',fontweight='bold')
    for ax, note in zip(axes,NOTES):
        for variant,label,color in zip(VARIANTS,LABELS,COLORS):
            angles = (0,45) if variant=='exterior' else (0,)
            for angle in angles:
                harmonics = rows[note,variant]['radiation'][str(angle)]['harmonics'][:6]
                values = [(h['order'],h['spl_db']) for h in harmonics if h['spl_db'] is not None]
                if values:
                    order,level = zip(*values)
                    suffix = f' ({angle}°)' if variant=='exterior' else ' (0° / 45°)'
                    ax.plot(order,level,marker='o',ms=4,lw=1.4,color=color,
                            ls='--' if angle else '-',label=label+suffix)
        ax.set(title=note.replace('Bb','B♭'),xlabel='Harmonic number',xticks=range(1,7))
        ax.grid(axis='y',alpha=.24)
    axes[0].set_ylabel('Harmonic level (dB re 20 µPa)')
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(handles,labels,fontsize=7,frameon=False)
    fig.savefig(out/'bell_harmonics.png',dpi=165)
    plt.close(fig)

    fig, (pitch_ax, cup_ax) = plt.subplots(1,2,figsize=(9,3.55),layout='constrained')
    x = np.arange(len(NOTES))
    for i,(variant,label,color) in enumerate(zip(VARIANTS,LABELS,COLORS)):
        pos = x + (i-1)*.12
        selected = [rows[note,variant] for note in NOTES]
        pitch_ax.plot(pos,[s['pitch_cents'] for s in selected],'o',color=color,label=label)
        cup_ax.plot(pos,[s['cup_rms_pa']/1000 for s in selected],'o',color=color)
    pitch_ax.axhline(0,lw=.7,color='#8b9299')
    pitch_ax.set(title='Sustained pitch',ylabel='Difference from target (cents)')
    cup_ax.set(title='Pressure inside mouthpiece cup',ylabel='Root mean square (kPa)')
    for ax in (pitch_ax,cup_ax):
        ax.set_xticks(x,[note.replace('Bb','B♭') for note in NOTES])
        ax.grid(axis='y',alpha=.24)
    pitch_ax.legend(fontsize=7,frameon=False)
    fig.savefig(out/'bell_pitch_and_pressure.png',dpi=165)
    plt.close(fig)


def write_wav(out, name, wave, gain):
    pcm = np.round(32767*gain*wave).astype('i2')
    wavfile.write(out/name,AUDIO_FS,pcm)


def write_audio(audio, out):
    """Publish common-gain clips so relative level and directivity stay audible."""
    gap = np.zeros(round(.15*AUDIO_FS))
    separation = np.zeros(round(.60*AUDIO_FS))
    peak = max(np.max(abs(angles[a])) for angles in audio.values() for a in (0,45))
    if not peak:
        raise ValueError('The comparison contains only silence.')
    gain = float(.89/peak)
    blocks = {}
    for variant in VARIANTS:
        for angle in ((0,45) if variant=='exterior' else (0,)):
            waves = [audio[note,variant][angle] for note in NOTES]
            blocks[variant,angle] = np.concatenate([waves[0],gap,waves[1],gap,waves[2]])
    order = []
    t = 0.
    for variant in VARIANTS:
        order.append(dict(variant=variant,start_s=t,angle_deg=0,note_order=list(NOTES)))
        t += (len(blocks[variant,0])+len(separation))/AUDIO_FS
    comparison = np.concatenate([blocks['previous',0],separation,blocks['fine_taper',0],
                                 separation,blocks['exterior',0]])
    write_wav(out,'bell_comparison_0deg.wav',comparison,gain)
    for angle in (0,45):
        write_wav(out,f'bell_exterior_{angle}deg.wav',blocks['exterior',angle],gain)
    return dict(sample_rate_hz=AUDIO_FS,common_gain_fs_per_pa=gain,comparison_order=order,
                comparison_duration_s=len(comparison)/AUDIO_FS,note_gap_s=len(gap)/AUDIO_FS)


def write_comparison(records, output_dir, *, metadata=None):
    """Export a complete three-note, three-variant comparison.

    Each record is keyed as (note, variant) and has the 192 kHz eleven-channel
    trace, raw 48 kHz audio in pascals keyed by angle, and target_hz. Optional
    keys include duration_s, the sustained window_s, and JSON diagnostics.
    """
    out = Path(output_dir)
    os.makedirs(out,exist_ok=True)
    rows, audio = {}, {}
    for note in NOTES:
        for variant in VARIANTS:
            key = note,variant
            rows[key],audio[key] = analyze_record(note,variant,records[key])
    result = dict(
        distance_m=1.,
        analysis='Audio and cup pitch estimated separately. Harmonic levels use a Hann window around the measured cup pitch.',
        notes=list(rows.values()), differences=comparisons(rows), audio=write_audio(audio,out),
        metadata=metadata or {},
    )
    plot_results(rows,out)
    (out/'bell_comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result
