"""Persist lossless states and full audit metadata for the new POM study."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import hashlib,json,os,zipfile
from scripts.run_pom_validation import ROOT,DOC,OUT,atomic_json


def main():
    registries=['community_registry.json','memory_registry.json','garden_registry.json','followup_registry.json','refinement_registry.json','rare_pool_registry.json']
    cases=[]
    for name in registries:cases+=json.loads((DOC/name).read_text())['cases']
    if len({c['id'] for c in cases})!=len(cases):raise ValueError('Duplicate archive IDs')
    files=[]
    for case in cases:
        for suffix in ['.npz','.json']:
            path=OUT/(case['id']+suffix)
            if not path.exists():raise ValueError('Missing '+str(path))
            files.append((path,'runs/'+path.name))
    for name in registries+['PROTOCOL.md','FOLLOWUP_PROTOCOL.md','REFINEMENT_PROTOCOL.md','RARE_POOL_PROTOCOL.md','invasion_diagnostic.json','EVIDENCE.md','RESULTS.md','summary.json','summary.csv','PATTERN_DECISIONS.json','VERIFICATION.json']:
        path=DOC/name
        if path.exists():files.append((path,'study/'+name))
    hashes=[dict(path=target,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path,target in files]
    contents=dict(study='SPLAS 0D POM challenge, 16 October 2026',cases=len(cases),
        initial_registered_community_cases=88,subsequent_density_controls=3,memory_cases=84,reserve_transplants=24,
        subsequent_whole_plane_and_corner_refinements=12,rare_reintroductions=10,
        units='Concentrations mmol element m^-3; days; reference size ESD micrometres. Abundance 1e8 cells m^-3.',
        state_format='Original exact float64 time/state/traits/grids, copied byte-for-byte; read-only flux observations reconstructed from states.',
        baseline_commit='8b14dda650cde76a9559e0d8ae73ae7bc3beb0a1',files=hashes)
    destination=ROOT/'output/data/SPLAS_Validazione_POM_0D_2026-10-16.zip';destination.parent.mkdir(parents=True,exist_ok=True)
    temp=destination.with_suffix('.zip.partial')
    readme='''SPLAS: independent 0D Pattern Oriented Modelling challenge\n\nThis archive contains only the NEW study. The earlier 117/166 scenarios are\nexploratory evidence and were not relabelled as confirmation. No parameter was\nfitted to the new pattern outcomes.\n\nRead study/RESULTS.md for successes, failures and scope; study/PROTOCOL.md was\nregistered internally before the primary runs. FOLLOWUP_PROTOCOL.md and REFINEMENT_PROTOCOL.md explicitly\nidentifies subsequent controls. The six registries contain exact parameters,\ninitial conditions, source SHA256 values and registration timestamps.\n\nEach runs/<ID>.npz has time, state, sizes, zsizes and traits. Its companion JSON\nrecords scientific sources, signature, archive SHA, budgets, steps, environment\nand exact common-state parents where applicable. NPZ bytes are unchanged.\nDiagnostics are cumulative flux ledgers, not live compartments. Source code and\nreproducible commands are in the repository's codex/dynamic-cn-development\nbranch and PR https://github.com/Marcofianchini/SPLAS-1.1/pull/1 . Check every\nscientific source against its registered SHA before reproducing a run.\n\nReproduction in that source version:\n  python scripts/run_pom_validation.py --workers 4\n  python scripts/run_pom_validation.py --stage memory --workers 4\n  python scripts/run_pom_validation.py --stage garden --workers 4\n  python scripts/run_pom_followup.py\n  python scripts/run_pom_refinement.py\n  python scripts/run_pom_rare_pool.py\n  python scripts/analyze_pom_invasion.py\n  python scripts/analyze_pom_validation.py\n  python scripts/plot_pom_validation.py\nNumba is required for the ensemble; Python/NumPy/SciPy/Numba versions are in\neach metadata record. Long runs use IMEX, and max_step controls splitting error.\nThe model and previous scientific sources are unmodified.\n\nThis is an ecological explanatory challenge. Conservation, numerical positivity\nand matching some qualitative patterns do not establish universal validation.\n'''
    with zipfile.ZipFile(temp,'w',allowZip64=True) as z:
        for path,target in files:
            z.write(path,target,compress_type=zipfile.ZIP_STORED if path.suffix=='.npz' else zipfile.ZIP_DEFLATED)
        z.writestr('CONTENTS.json',json.dumps(contents,indent=2)+'\n',compress_type=zipfile.ZIP_DEFLATED)
        z.writestr('README.txt',readme,compress_type=zipfile.ZIP_DEFLATED)
    with temp.open('rb') as f:os.fsync(f.fileno())
    with zipfile.ZipFile(temp) as z:
        bad=z.testzip()
        if bad:raise IOError('ZIP verification failed: '+bad)
    temp.replace(destination)
    record=dict(**contents,archive_path=str(destination.relative_to(ROOT)),bytes=destination.stat().st_size,
        sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),zip_crc_verified=True)
    atomic_json(DOC/'data_manifest.json',record)
    print(json.dumps({k:record[k] for k in ['cases','archive_path','bytes','sha256','zip_crc_verified']}))

if __name__=='__main__':main()
