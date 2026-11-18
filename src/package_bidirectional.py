"""Package reproducible sources, final presets, audio, and measured diagnostics."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import zipfile

ROOT=Path(__file__).resolve().parents[1]


def main(destination=None):
    folder=ROOT/'output/bidirectional'
    report=json.loads((folder/'verification.json').read_text())
    assert report['passed'] and json.loads((folder/'solver_tests.json').read_text())['passed']
    files=[ROOT/name for name in ['README.md','docs/ja/NONLINEAR_MODEL.md',
                                 'docs/ja/BIDIRECTIONAL_MODEL.md','requirements.txt']]
    files.extend(ROOT/'src'/name for name in ['trumpet.py','nonlinear_propagation.py','bidirectional.py',
        'render_bidirectional.py','verify_bidirectional.py','plot_bidirectional.py','package_bidirectional.py'])
    files.extend(ROOT/'data'/name for name in ['model.json','performance.json','performance_bidirectional.json',
        'york_trumpet_fig7_52.csv','sources.json','sources_bidirectional.json'])
    files.extend(folder/name for name in ['trumpet_Bb3_Bb5_major_f_bidirectional_nonlinear.wav',
        'trumpet_Bb3_Bb5_major_f_bidirectional_linear_reference.wav',
        'comparison_bidirectional_linear_then_nonlinear.wav','measurements.json',
        'verification.json','solver_tests.json','diagnostics.png'])
    manifest={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    archive=ROOT/'trumpet_Bb3_Bb5_bidirectional_physical_model.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files:z.write(p,str(p.relative_to(ROOT)))
        z.writestr('SHA256.json',json.dumps(manifest,indent=2))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for name,digest in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==digest
    if destination:
        dest=Path(destination);os.makedirs(dest,exist_ok=True)
        exports=[archive,ROOT/'docs/ja/BIDIRECTIONAL_MODEL.md',folder/'diagnostics.png',folder/'verification.json',
            folder/'trumpet_Bb3_Bb5_major_f_bidirectional_nonlinear.wav',
            folder/'trumpet_Bb3_Bb5_major_f_bidirectional_linear_reference.wav',
            folder/'comparison_bidirectional_linear_then_nonlinear.wav']
        for p in exports:
            target=dest/p.name
            # Name generic diagnostics explicitly to preserve previous versions.
            if p.name in ['diagnostics.png','verification.json']:target=dest/('bidirectional_'+p.name)
            shutil.copy2(p,target)
            assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256(target.read_bytes()).digest()
    print(json.dumps(dict(archive=str(archive),files=len(files)+1,bytes=archive.stat().st_size,destination=destination),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--dest');args=parser.parse_args();main(args.dest)
