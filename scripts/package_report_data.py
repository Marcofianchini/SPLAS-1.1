"""Attach generated numerical evidence to the PDF; source code stays in Git.

Run only after final strict analysis and a fresh PDF render. The attachment
contains ignored generated results, never repository source or user thesis PDFs.
"""
from pathlib import Path
import io,json,hashlib,zipfile,argparse
from pypdf import PdfReader,PdfWriter
ROOT=Path(__file__).resolve().parents[1]
PDF=ROOT/'output/pdf/SPLAS_Sviluppo_Verifica_2026-10-16.pdf'
NAME='SPLAS_results_2026-10-16.zip'

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--allow-partial',action='store_true',help='Preserve an explicitly incomplete report checkpoint');args=ap.parse_args()
 summary=json.load(open(ROOT/'research/validation_summary.json'))
 missing=set(summary['missing_cases'])
 if missing and not args.allow_partial:raise RuntimeError('Incomplete numerical evidence')
 registry=[]
 for name in ['experiment_registry.json','supplement_registry.json','carbon_boundary_registry.json','recommended_registry.json','recommended_refinement_registry.json','export_intervention_registry.json']:
  p=ROOT/'research'/name
  if p.exists():registry+=json.load(open(p))
 files=[]
 for case in registry:
  if case['id'] in missing:continue
  for suffix in ['.json','.npz']:files.append(ROOT/'results/thesis_replay'/(case['id']+suffix))
 files+=sorted((ROOT/'results/thesis_replay').glob('run_status*.json'))
 for directory in ['cultures','mechanisms','protocols','figures','recommended_initials']:
  files+=sorted((ROOT/'results'/directory).glob('*'))
 files=[p for p in files if p.is_file() and p.suffix in ['.json','.npz','.png','.pdf']]
 manifest={'source_repository':'https://github.com/Marcofianchini/SPLAS-1.1','branch':'codex/dynamic-cn-development',
  'model_sha256':summary['model_sha256'],'data_citation':'Liefer etal2026,PANGAEA900248,CC-BY4.0; article10.3389/fmicb.2019.00763',
  'interpretation':'Generated numerical evidence only; no field/export validation. Source code and tracked configuration files are versioned separately in Git.',
  'complete':not bool(missing),'missing_cases':sorted(missing),'files':[]}
 stream=io.BytesIO()
 with zipfile.ZipFile(stream,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
  for p in files:
   data=p.read_bytes();name=p.relative_to(ROOT/'results').as_posix();archive.writestr(name,data)
   manifest['files'].append({'path':name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
  archive.writestr('manifest.json',json.dumps(manifest,indent=2))
  archive.writestr('READ_ME.txt','SPLAS-CN generated numerical evidence.\nAllC:N ratios are molar; chemical unitsmmol element/m3; time days.\nModel source and assumptions:GitHub branch codex/dynamic-cn-development.\nFiles are generated results, not experimental observations except culture observation arrays marked observed.\nSee the report for failed predictions and validation limits.\n')
 payload=stream.getvalue();reader=PdfReader(PDF)
 if reader.attachments:raise RuntimeError('Regenerate an unattached base report before packaging')
 writer=PdfWriter();writer.clone_document_from_reader(reader);writer.add_attachment(NAME,payload)
 writer.add_metadata({'/Subject':'Research report with embedded numerical results archive; C:N verification and exploratory validation.'})
 temp=PDF.with_suffix('.tmp.pdf')
 with temp.open('wb') as out:writer.write(out)
 temp.replace(PDF)
 check=PdfReader(PDF);attachments=check.attachments
 assert NAME in attachments and len(attachments[NAME])==1 and hashlib.sha256(attachments[NAME][0]).digest()==hashlib.sha256(payload).digest()
 record={'attachment':NAME,'attachment_bytes':len(payload),'attachment_sha256':hashlib.sha256(payload).hexdigest(),'generated_files':len(files),'pdf_pages':len(check.pages),'pdf_bytes':PDF.stat().st_size,'pdf_sha256':hashlib.sha256(PDF.read_bytes()).hexdigest()}
 (ROOT/'research/report_package_manifest.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
if __name__=='__main__':main()
