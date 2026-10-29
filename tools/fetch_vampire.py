"""Fetch the exact public VAMPIRE base into an empty destination."""
from pathlib import Path
import argparse,json,subprocess
ROOT=Path(__file__).resolve().parents[1]
def fetch(destination):
    destination=Path(destination).resolve();spec=json.loads((ROOT/'upstream.json').read_text())['vampire']
    if destination.exists():raise SystemExit('Destination exists; use a new directory')
    destination.parent.mkdir(parents=True,exist_ok=True)
    subprocess.run(['git','-c','core.autocrlf=false','clone','--depth','1','--branch',spec['tag'],spec['repository']+'.git',str(destination)],check=True)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=destination,text=True).strip()
    if head!=spec['commit']:raise SystemExit('Upstream commit mismatch')
    print('Verified VAMPIRE base:',head)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--destination',type=Path,default=ROOT/'build/vampire');a=p.parse_args();fetch(a.destination)
