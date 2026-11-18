"""Tune lip spring frequency, never resample the produced pitch."""
from bidirectional import *


def run():
    settings=json.loads((ROOT/'data/performance.json').read_text())
    result=[];log=[];cache={}
    for s in settings:
        valves=tuple(s['valves_semitones'])
        if valves not in cache:cache[valves]=setup(valves)
        cfg=cache[valves];target=s['target_hz'];initial=s['lip_frequency_hz']
        duration=1.0 if s['note']=='Bb5' else .65
        trials=[]
        def evaluate(flip):
            tr,mach=simulate_coupled(cfg,duration,flip,s['blowing_pressure_pa'],s['rest_y_m'],hold=duration)
            m=measure(tr,fs=cfg['fs'],start=duration-.22)
            cents=float(1200*np.log2(m['frequency']/target))
            trials.append(dict(lip_frequency_hz=flip,cents=cents,**m))
            print(s['note'],round(flip,5),round(m['frequency'],5),round(cents,3),flush=True)
            return m['frequency']
        x0=initial;y0=evaluate(x0)
        best=(abs(1200*np.log2(y0/target)),x0)
        if best[0]>1.0:
            x1=x0+max(2.,initial*.007);y1=evaluate(x1)
            best=min(best,(abs(1200*np.log2(y1/target)),x1))
            for _ in range(5):
                if best[0]<.7:break
                slope=(y1-y0)/(x1-x0)
                if slope<=.02:break
                step=np.clip((target-y1)/slope,-15.,15.)
                x2=float(np.clip(x1+step,initial*.90,initial*1.12))
                if abs(x2-x1)<.01:break
                y2=evaluate(x2)
                best=min(best,(abs(1200*np.log2(y2/target)),x2))
                x0,y0,x1,y1=x1,y1,x2,y2
        s['lip_frequency_hz']=best[1]
        result.append(s);log.append(dict(note=s['note'],trials=trials,best_cents=best[0]))
        (ROOT/'data/performance_bidirectional.json').write_text(json.dumps(result,indent=2))
        (ROOT/'output/bidirectional/tuning.json').write_text(json.dumps(log,indent=2))


if __name__=='__main__':run()
