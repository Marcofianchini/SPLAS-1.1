"""Analytical closure scaling control, declared after the primary plane."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import json,hashlib
from concurrent.futures import ProcessPoolExecutor,as_completed
from scripts.run_pom_validation import registry,identifier,execute,atomic_json,DOC,source_hashes,BASELINE
from datetime import datetime,timezone


def cases():
    by={c['id']:c for c in registry()};result=[]
    for light,n in [(3.,.5),(300.,.5),(300.,8.)]:
        c=json.loads(json.dumps(by[identifier('own',light,n)+'_bins70']))
        c.update(id=c['id']+'_density',family='closure_density_control',
            parent_reference=identifier('own',light,n),normalization_factor=69/34,
            followup_protocol_sha256=hashlib.sha256((DOC/'FOLLOWUP_PROTOCOL.md').read_bytes()).hexdigest())
        c['physiology']['zoo_mortality']*=69/34;result.append(c)
    return result


def main():
    cc=cases();path=DOC/'followup_registry.json'
    reg=dict(baseline_commit=BASELINE,registered_before_integration=True,
        followup_protocol_sha256=hashlib.sha256((DOC/'FOLLOWUP_PROTOCOL.md').read_bytes()).hexdigest(),
        source_hashes=source_hashes(),cases=cc)
    if path.exists():
        old=json.loads(path.read_text());old.pop('registered_at_UTC',None)
        if reg!=old:raise ValueError('Changed followup registry')
    else:
        reg['registered_at_UTC']=datetime.now(timezone.utc).isoformat();atomic_json(path,reg)
    print('Registered 3 closure-density controls before their integration',flush=True)
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(execute,c):c['id'] for c in cc}
        for future in as_completed(futures):
            r=future.result();print(json.dumps(dict(id=r['case']['id'],seconds=r['runtime_seconds'],residual=r['maximum_relative_residual'])),flush=True)

if __name__=='__main__':main()
