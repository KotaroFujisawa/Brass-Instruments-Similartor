"""Find self-oscillating embouchures and render a two-octave Bb-major scale."""
from trumpet import *
from scipy.io import wavfile
import argparse
import time
import hashlib

NOTES=[('Bb3',58,()),('C4',60,(2,3)),('D4',62,(2,1)),('Eb4',63,(2,)),
       ('F4',65,()),('G4',67,(2,1)),('A4',69,(1,)),('Bb4',70,()),
       ('C5',72,(2,)),('D5',74,()),('Eb5',75,(2,)),('F5',77,()),
       ('G5',79,(2,1)),('A5',81,(1,)),('Bb5',82,())]


def evaluate(res,flip,pressure,yeq,duration=.48):
    tr=simulate(res['r'],float(res['z0']),FS,duration,flip,pressure,yeq,hold=duration)
    m=measure(tr,start=duration-.17)
    m.update(lip_frequency_hz=float(flip),blowing_pressure_pa=float(pressure),rest_y_m=float(yeq))
    return m


def fit_note(name,target,res):
    pressure=5000*(target/233.08188)**.60
    ys=[0.,-.00004,-.00008,.00002] if target>620 else [0.,.00002,-.00004]
    trials=[]
    def attempt(f,p,y):
        try:
            m=evaluate(res,f,p,y)
            if m['frequency']>0:
                m['cents']=float(1200*np.log2(m['frequency']/target))
                trials.append(m)
                return m
        except ValueError: pass
        return None
    # Search resonance-capture windows; genuinely self-excited, with no forced tone.
    for y in ys:
        scan=[]
        for ratio in np.linspace(.96,1.26,9):
            m=attempt(target*ratio,pressure,y)
            if m is not None and m['rms']>400 and m['periodicity']>.98:
                scan.append(m)
        for a,b in zip(scan[:-1],scan[1:]):
            if a['cents']*b['cents']<=0 and abs(a['cents']-b['cents'])<180:
                lo,hi=a,b
                for _ in range(7):
                    mid=attempt((lo['lip_frequency_hz']+hi['lip_frequency_hz'])/2,pressure,y)
                    if mid is None or mid['rms']<300 or abs(mid['cents'])>150:break
                    if abs(mid['cents'])<.3:break
                    if mid['cents']*lo['cents']<=0:hi=mid
                    else:lo=mid
        valid=[m for m in trials if m['rms']>500 and m['periodicity']>.985 and abs(m['cents'])<2]
        if valid:
            best=min(valid,key=lambda m:abs(m['cents']))
            print(name,'FIT',best,flush=True)
            return best,trials
    # Increased pressure within an explicitly assumed forte/high-register range.
    for factor in [1.2,1.4]:
        p=pressure*factor
        for y in [-.00004,0.]:
            scan=[]
            for ratio in np.linspace(1.0,1.27,9):
                m=attempt(target*ratio,p,y)
                if m and m['rms']>500 and m['periodicity']>.985:scan.append(m)
            for a,b in zip(scan[:-1],scan[1:]):
                if a['cents']*b['cents']<=0 and abs(a['cents']-b['cents'])<150:
                    lo,hi=a,b
                    for _ in range(8):
                        mid=attempt((lo['lip_frequency_hz']+hi['lip_frequency_hz'])/2,p,y)
                        if mid is None or mid['rms']<300 or abs(mid['cents'])>150:break
                        if abs(mid['cents'])<.3:break
                        if mid['cents']*lo['cents']<=0:hi=mid
                        else:lo=mid
            good=[m for m in trials if m['rms']>500 and m['periodicity']>.985 and abs(m['cents'])<2]
            if good:
                best=min(good,key=lambda m:abs(m['cents']))
                print(name,'FIT',best,flush=True)
                return best,trials
    good=[m for m in trials if m['rms']>500 and m['periodicity']>.98]
    if not good:raise RuntimeError(f'No stable oscillation: {name}')
    best=min(good,key=lambda m:abs(m['cents']))
    print(name,'NEEDS_REVIEW',best,flush=True)
    return best,trials


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fit',action='store_true')
    parser.add_argument('--render',action='store_true')
    args=parser.parse_args()
    (ROOT/'tmp').mkdir(exist_ok=True)
    (ROOT/'output').mkdir(exist_ok=True)
    params=json.loads((ROOT/'data/model.json').read_text())
    dep=params['cup_depth_m']
    signature=hashlib.sha256((ROOT/'data/model.json').read_bytes()
        +(ROOT/'data/york_trumpet_fig7_52.csv').read_bytes()
        +(ROOT/'src/trumpet.py').read_bytes()).hexdigest()
    cache={}
    for _,_,valves in NOTES:
        if valves in cache:continue
        file=ROOT/'tmp'/('res_'+''.join(map(str,valves))+'.npz')
        res=dict(np.load(file)) if file.exists() else {}
        if str(res.get('signature',''))!=signature:
            res=resonator(valves,dep)
            res['signature']=signature
            np.savez(file,**res)
        cache[valves]=res
        print('Resonator',valves,'peaks',res['peaks'][:10],flush=True)
    settings=[]
    if args.fit:
        for name,midi,valves in NOTES:
            dest=ROOT/'tmp'/f'fit_{name}.json'
            target=440*2**((midi-69)/12)
            saved=json.loads(dest.read_text()) if dest.exists() else {}
            if saved.get('signature')==signature:
                best=saved['selected']
            else:
                best,trials=fit_note(name,target,cache[valves])
                dest.write_text(json.dumps(dict(selected=best,trials=trials,signature=signature),indent=2))
            best.update(note=name,midi=midi,target_hz=target,valves_semitones=list(valves))
            settings.append(best)
        (ROOT/'data/performance.json').write_text(json.dumps(settings,indent=2))
    if args.render:
        settings=json.loads((ROOT/'data/performance.json').read_text())
        waves=[];offset=0.;measurements=[]
        for i,s in enumerate(settings):
            res=cache[tuple(s['valves_semitones'])]
            duration=1.1 if i==len(settings)-1 else .72
            tr=simulate(res['r'],float(res['z0']),FS,duration,s['lip_frequency_hz'],
                        s['blowing_pressure_pa'],s['rest_y_m'],hold=duration-.14)
            m=measure(tr,start=.28,end=duration-.15)
            m['cents']=float(1200*np.log2(m['frequency']/s['target_hz'])) if m['frequency']>0 else None
            m.update(note=s['note'],start_s=offset,duration_s=duration,target_hz=s['target_hz'])
            measurements.append(m)
            y=radiate(tr,res)
            # Antialiasing filter is the only tonal processing of the radiated pressure.
            y=signal.resample_poly(y,1,2)
            waves.append(y);waves.append(np.zeros(int(.06*48000)))
            offset+=duration+.06
            np.savez_compressed(ROOT/'output'/f'trace_{s["note"]}.npz',trace=tr[::2],fs=48000,
                                radiation_pa=y,settings=json.dumps(s))
            print('Render',m,flush=True)
        audio=np.concatenate(waves)
        # One gain for the whole scale, not independent loudness normalisation per note.
        gain=.89/np.max(abs(audio))
        pcm=np.round(audio*gain*32767).astype(np.int16)
        wavfile.write(ROOT/'output/trumpet_Bb3_Bb5_major_f.wav',48000,pcm)
        wavfile.write(ROOT/'output/trumpet_Bb3_Bb5_radiation_Pa_float.wav',48000,audio.astype(np.float32))
        (ROOT/'output/measurements.json').write_text(json.dumps(dict(
            sample_rate=48000,duration_s=len(audio)/48000,common_gain_fs_per_pa=gain,
            peak_radiated_pa=float(np.max(abs(audio))),notes=measurements),indent=2))
        print('WAV completed:',len(audio)/48000,'seconds',flush=True)


if __name__=='__main__':main()
