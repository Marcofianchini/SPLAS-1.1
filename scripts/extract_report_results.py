"""Extract the numerical ZIP embedded in the distributed validation PDF."""
from pathlib import Path,PurePosixPath
import argparse,io,zipfile,json,hashlib
from pypdf import PdfReader

def main():
 ap=argparse.ArgumentParser();ap.add_argument('pdf',type=Path);ap.add_argument('--output',type=Path,default=Path('report_results'));args=ap.parse_args()
 reader=PdfReader(args.pdf);name='SPLAS_results_2026-10-16.zip';items=reader.attachments.get(name)
 if not items or len(items)!=1:raise RuntimeError('Expected one SPLAS numerical archive')
 with zipfile.ZipFile(io.BytesIO(items[0])) as z:
  for n in z.namelist():
   p=PurePosixPath(n)
   if p.is_absolute() or '..' in p.parts:raise ValueError('Unsafe archive path')
  args.output.mkdir(parents=True,exist_ok=True);z.extractall(args.output)
 manifest=json.load(open(args.output/'manifest.json'))
 for item in manifest['files']:
  content=(args.output/item['path']).read_bytes()
  if len(content)!=item['bytes'] or hashlib.sha256(content).hexdigest()!=item['sha256']:raise RuntimeError('Checksum mismatch: '+item['path'])
 print(f"Extracted and verified {len(manifest['files'])} generated files in {args.output.resolve()}")
if __name__=='__main__':main()
