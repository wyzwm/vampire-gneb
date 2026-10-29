"""Build only the serial extension executable from a verified installed tree."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
def build(source,jobs,compiler,make):
    source=Path(source).resolve();marker=source/'.gneb-extension.json'
    if not marker.exists():raise SystemExit('Run install_extension.py first')
    record=json.loads(marker.read_text())
    for name,h in record['installed_files'].items():
        if hashlib.sha256((source/name).read_bytes()).hexdigest()!=h:raise SystemExit('Installed source changed: '+name)
    if not 1<=jobs<=64:raise SystemExit('Choose --jobs between 1 and 64')
    out=ROOT/'build/bin';out.mkdir(parents=True,exist_ok=True)
    exe=out/('vampire-gneb.exe' if os.name=='nt' else 'vampire-gneb')
    cxx=shutil.which(compiler) if not Path(compiler).is_file() else str(Path(compiler).resolve())
    maker=shutil.which(make) if not Path(make).is_file() else str(Path(make).resolve())
    if not cxx or not maker:raise SystemExit('C++ compiler or make not found; see docs/INSTALL.md')
    env=os.environ.copy();env['PATH']=str(Path(cxx).parent)+os.pathsep+str(Path(maker).parent)+os.pathsep+env.get('PATH','')
    # Use only the basename in make's shell command, retaining the configured executable on PATH.
    cmd=[maker,'serial','-j'+str(jobs),'GCC='+Path(cxx).name+' -std=c++11'+(' -DWIN_COMPILE' if os.name=='nt' else ''),
         'GCC_CFLAGS=-O2 -std=c++11 -D_USE_MATH_DEFINES -I./hdr -I./src/qvoronoi',
         'GCC_LDFLAGS='+('-static -static-libgcc -static-libstdc++' if os.name=='nt' else '-lstdc++'),
         'GHASH='+record['upstream_commit']+'+gneb-'+record['extension_version'],
         'EXECUTABLE='+os.path.relpath(exe,source).replace('\\','/')]
    # Upstream tracks empty object directories; support source archives too.
    for directory in (source/'src').rglob('*'):
        if directory.is_dir():(source/'obj'/directory.relative_to(source/'src')).mkdir(parents=True,exist_ok=True)
    log=ROOT/'build/build.log';start=time.monotonic()
    with log.open('wb') as f:r=subprocess.run(cmd,cwd=source,env=env,stdout=f,stderr=subprocess.STDOUT)
    result={'exit_code':r.returncode,'elapsed_s':time.monotonic()-start,'compiler':subprocess.check_output([cxx,'--version'],env=env,text=True).splitlines()[0],
            'sha256':hashlib.sha256(exe.read_bytes()).hexdigest() if r.returncode==0 else None,'serial_only':True}
    (ROOT/'build/build-result.json').write_text(json.dumps(result,indent=2)+'\n')
    if r.returncode:
        print('\n'.join(log.read_text(errors='replace').splitlines()[-25:]));raise SystemExit(r.returncode)
    print('Built',exe);print('SHA256',result['sha256'])
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'build/vampire');p.add_argument('--jobs',type=int,default=2);p.add_argument('--compiler',default='g++');p.add_argument('--make',default='make');a=p.parse_args();build(a.source,a.jobs,a.compiler,a.make)
