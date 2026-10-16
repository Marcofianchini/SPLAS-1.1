"""Lossless state archive; regenerate read-only observations without rerunning ODEs.

The archive contains every saved state in float64 (including cumulative flux
ledgers), time, traits, grids, original run metadata and hashes. Bulky process
arrays are deterministic observations of these states and stay reproducible
with the matching Git revision. No repository or third-party source is bundled.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse,hashlib,io,json,zipfile
from dataclasses import replace
import numpy as np
from splas2 import Model
from splas2.stoichiometry import observe
from scripts.run_stoichiometry import ROOT,OUT,registry,build,signature
from scripts.run_stoichiometry_controls import registry as controls
from scripts.run_stoichiometry_robustness import registry as robustness

KEYS=['time','state','sizes','zsizes','traits']

def digest(data):return hashlib.sha256(data).hexdigest()

def restore_arrays(case,base):
    m=build(case);observations={}
    for t,y in zip(base['time'],base['state']):
        mm=m
        if case['segments']:
            boundary=0.
            for duration,light,pulse in case['segments']:
                if t<boundary+duration-1e-9 or boundary+duration>=case['days']:
                    mm=Model(replace(m.config,surface_PAR=light),m.physiology);break
                boundary+=duration
        for key,value in observe(mm,y,float(t)).items():
            observations.setdefault(key,[]).append(value)
    return {**base,**{k:np.asarray(v) for k,v in observations.items()}}

def pack(destination):
    destination.parent.mkdir(parents=True,exist_ok=True);cases=registry()+controls()+robustness();manifest=[]
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_STORED) as archive:
        for case in cases:
            id=case['id'];path=OUT/(id+'.npz');metadata=json.loads(path.with_suffix('.json').read_text())
            if metadata['run_signature']!=signature(case) or metadata['archive_sha256']!=digest(path.read_bytes()):
                raise ValueError('Unverified source '+id)
            with np.load(path) as full:
                base={k:full[k] for k in KEYS}
                # Check representative temporal, uniform and resolved communities
                # bitwise against all originally saved derived observations.
                if id in ['T_N_then_light','C_L45_N1_uniform','R_L405_N1_bins140']:
                    rebuilt=restore_arrays(case,base)
                    for key,value in rebuilt.items():
                        if not np.array_equal(value,full[key]):raise ValueError('Nonidentical reconstructed array '+id+':'+key)
            buf=io.BytesIO();np.savez_compressed(buf,**base);data=buf.getvalue()
            archive.writestr('states/'+id+'.npz',data)
            archive.writestr('metadata/'+id+'.json',json.dumps(metadata,indent=2)+'\n')
            manifest.append(dict(id=id,compact_sha256=digest(data),original_archive_sha256=metadata['archive_sha256'],
                arrays={k:dict(shape=list(v.shape),dtype=str(v.dtype),sha256=digest(v.tobytes())) for k,v in base.items()}))
        for name in ['registry.json','control_registry.json','robustness_registry.json','summary.json','summary.csv','DESIGN.md','RESULTS.md','LITERATURE.md','PATTERNS.json']:
            source=ROOT/'research/stoichiometry'/name
            if source.exists():archive.write(source,'study/'+name)
        archive.writestr('MANIFEST.json',json.dumps(dict(schema=1,lossless_states=True,cases=manifest),indent=2)+'\n')
        archive.writestr('README.txt','SPLAS: emergent ecological patterns from dynamic C/N\n'
            'All saved states are exact float64, with integrated C/N process ledgers.\n'
            'Derived instantaneous process arrays were omitted to avoid redundant storage.\n'
            'Matching code: github.com/Marcofianchini/SPLAS-1.1 branch codex/dynamic-cn-development.\n'
            'Source hashes and exact environments are in each run metadata file.\n'
            'Restore: python scripts/package_stoichiometry.py --restore ARCHIVE.zip\n'
            'Analyze: python scripts/analyze_stoichiometry.py\n'
            'No ODE integration is needed to reconstruct the process arrays.\n')
    metadata=dict(file_name=destination.name,bytes=destination.stat().st_size,sha256=digest(destination.read_bytes()),
        cases=len(cases),stored_arrays=KEYS,omitted='Instantaneous processes, exactly recomputable from states and matching code',
        reconstruction_checks=['T_N_then_light','C_L45_N1_uniform','R_L405_N1_bins140'])
    (ROOT/'research/stoichiometry/data_manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata,indent=2))

def restore(source):
    OUT.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(source) as archive:
        manifest=json.loads(archive.read('MANIFEST.json'))
        for record in manifest['cases']:
            id=record['id'];path=OUT/(id+'.npz');jp=path.with_suffix('.json')
            if path.exists() or jp.exists():raise FileExistsError('Refusing to replace existing run '+id)
            data=archive.read('states/'+id+'.npz')
            if digest(data)!=record['compact_sha256']:raise ValueError('Corrupt compact state '+id)
            with np.load(io.BytesIO(data)) as z:base={k:z[k] for k in KEYS}
            for k,check in record['arrays'].items():
                if digest(base[k].tobytes())!=check['sha256']:raise ValueError('Corrupt array '+id+':'+k)
            metadata=json.loads(archive.read('metadata/'+id+'.json'));case=metadata['case']
            if metadata['run_signature']!=signature(case):raise ValueError('Use matching model/observation source hashes: '+id)
            np.savez_compressed(path,**restore_arrays(case,base))
            metadata['restored_from_lossless_state_archive']=str(Path(source).name)
            metadata['original_archive_sha256']=metadata['archive_sha256']
            metadata['archive_sha256']=digest(path.read_bytes());jp.write_text(json.dumps(metadata,indent=2)+'\n')
            print('Restored '+id,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--restore',type=Path);a=parser.parse_args()
    if a.restore:restore(a.restore)
    else:pack(ROOT/'output/data/SPLAS_Patterns_Stechiometria_2026-10-16.zip')
