"""Small reproducible entry point; full thesis ensembles live in scripts/."""
import argparse,json
from pathlib import Path
from .model import *

def main():
 ap=argparse.ArgumentParser(description='Conservative size-structured plankton C:N model')
 ap.add_argument('--config',type=Path,help='JSON containing config, physiology and optional initial/max_step/method/days')
 ap.add_argument('--days',type=float);ap.add_argument('--output',type=Path,default=Path('results/example.npz'))
 ap.add_argument('--method',choices=['explicit','imex']);args=ap.parse_args()
 spec=json.loads(args.config.read_text()) if args.config else {}
 m=Model(Config(**spec.get('config',{'n_phyto':20,'n_zoo':20})),Physiology(**spec.get('physiology',{})))
 days=args.days if args.days is not None else spec.get('days',30.)
 y=m.initial(**spec.get('initial',{}))
 r=simulate(m,days,y0=y,max_step=spec.get('max_step',.05),method=args.method or spec.get('method','imex'),save_interval=min(10.,days/50.))
 r.save(args.output);print(json.dumps(r.diagnostics(),indent=2));print(args.output.resolve())
if __name__=='__main__':main()
