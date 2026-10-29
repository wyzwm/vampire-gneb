"""Install the reviewed overlay into a clean, exact VAMPIRE 7.0 checkout.
No compilation, simulation, network request or production executable replacement.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,shutil
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def normalized_sha(p):return hashlib.sha256(Path(p).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def install(source):
    source=Path(source).resolve();manifest=json.loads((ROOT/'upstream.json').read_text());marker=source/'.gneb-extension.json'
    if marker.exists():
        record=json.loads(marker.read_text())
        if record['patch_sha256']!=manifest['patch_sha256']:raise SystemExit('A different extension version is installed')
        for path,digest in record['installed_files'].items():
            if sha(source/path)!=digest:raise SystemExit('Installed file changed: '+path)
        print('Exact extension already installed; nothing changed.');return
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip()
    if head!=manifest['vampire']['commit']:raise SystemExit('Requires exact VAMPIRE 7.0.0 commit '+manifest['vampire']['commit'])
    dirty=subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=source,text=True)
    if dirty.strip():raise SystemExit('Tracked source changes exist; use a clean separate checkout')
    for name,digest in manifest['upstream_base_files'].items():
        if normalized_sha(source/name)!=digest:raise SystemExit('Base-file mismatch: '+name)
    for name,digest in manifest['overlay_sha256'].items():
        if sha(ROOT/'extension'/name)!=digest:raise SystemExit('Overlay checksum mismatch: '+name)
        if (source/name).exists():raise SystemExit('Overlay target already exists: '+name)
    patch=ROOT/'patches/vampire-7.0.0-integration.patch'
    if sha(patch)!=manifest['patch_sha256']:raise SystemExit('Patch checksum mismatch')
    subprocess.run(['git','apply','--check',str(patch)],cwd=source,check=True)
    subprocess.run(['git','apply',str(patch)],cwd=source,check=True)
    paths=list(manifest['upstream_base_files'])
    for name in manifest['overlay_sha256']:
        target=source/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/'extension'/name,target);paths.append(name)
    for origin,name in [('NOTICES.md','GNEB_NOTICES.md'),('LICENSES/Spirit-MIT.txt','SPIRIT_LICENSE.txt')]:
        shutil.copyfile(ROOT/origin,source/name);paths.append(name)
    record={'extension_version':manifest['extension_version'],'upstream_commit':head,'patch_sha256':manifest['patch_sha256'],
            'installed_files':{name:sha(source/name) for name in paths}}
    with marker.open('x',encoding='utf-8') as f:json.dump(record,f,indent=2);f.write('\n')
    print('Installed serial GNEB extension; source hashes recorded in .gneb-extension.json')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);a=p.parse_args();install(a.source)
