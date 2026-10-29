"""Public-data-only acceptance suite. Standard library plus a C++11 compiler."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,hashlib,importlib.util,json,math,os,shutil,subprocess,tempfile
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('example_runner',ROOT/'tools/run_example.py');example=importlib.util.module_from_spec(spec);spec.loader.exec_module(example)
def spins(path):
    lines=path.read_text().splitlines();return [list(map(float,l.split())) for l in lines[2:]]
def write_spins(path,rows):
    with path.open('w',encoding='ascii') as f:f.write('GNEB_SPINS_V1\n'+str(len(rows))+'\n');f.writelines(' '.join(format(x,'.17g') for x in r)+'\n' for r in rows)
def main(exe,cxx,report,threads=1):
    exe=Path(exe).resolve();compiler=shutil.which(cxx) or str(Path(cxx).resolve())
    (ROOT/'build').mkdir(exist_ok=True);out=Path(tempfile.mkdtemp(prefix='acceptance-',dir=ROOT/'build'))
    core=out/('test-core.exe' if os.name=='nt' else 'test-core')
    subprocess.run([compiler,'-std=c++11','-O2','-Wall','-Wextra','-pedantic']+(['-static'] if os.name=='nt' else [])+[str(ROOT/'tests/test_core.cpp'),'-o',str(core)],check=True)
    coretest=subprocess.run([str(core)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60);assert coretest.returncode==0,coretest.stderr.decode(errors='replace')
    result9=example.run(exe,out/'public9')
    source=ROOT/'examples/spirit-skyrmion-16x16'
    def copy_case(name):
        d=out/name;d.mkdir()
        for f in source.iterdir():
            if f.name in ['input','Co.mat','file.ucf','chain.txt','provenance.json'] or f.name.startswith('image_'):shutil.copyfile(f,d/f.name)
        return d
    d17=copy_case('public17');first=spins(d17/'image_0.txt');last=spins(d17/'image_8.txt')
    for k in range(17):
        fraction=k/16;rows=[]
        for a,b in zip(first,last):
            sa=a[3:];sb=b[3:];dot=max(-1.,min(1.,sum(x*y for x,y in zip(sa,sb))));theta=math.acos(dot)
            if theta<1e-9:s=[(1-fraction)*x+fraction*y for x,y in zip(sa,sb)]
            else:s=[(math.sin((1-fraction)*theta)*x+math.sin(fraction*theta)*y)/math.sin(theta) for x,y in zip(sa,sb)]
            norm=math.sqrt(sum(x*x for x in s));rows.append(a[:3]+[x/norm for x in s])
        write_spins(d17/f'resampled_{k}.txt',rows)
    (d17/'chain.txt').write_text(''.join(f'resampled_{k}.txt\n' for k in range(17)))
    r=subprocess.run([str(exe)],cwd=d17,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180);assert r.returncode==0,r.stderr.decode(errors='replace')
    result17=example.verify(d17);assert abs(result17['forward_barrier_meV']-result9['forward_barrier_meV'])<1e-5
    negative=[]
    for name,expected in [('nonzero_temperature',1),('unrelaxed_endpoint',1),('duplicate_image',1),('antipodal_neighbour',1),('wrong_coordinates',1),('nonunit_spin',1),('iteration_budget',2)]:
        d=copy_case(name)
        if name=='nonzero_temperature':(d/'input').write_text((d/'input').read_text().replace('sim:temperature=0','sim:temperature=2'))
        elif name=='unrelaxed_endpoint':
            rows=spins(d/'image_0.txt');rows[0][3:]=[1.,0.,0.];write_spins(d/'image_0.txt',rows)
        elif name=='duplicate_image':shutil.copyfile(d/'image_2.txt',d/'image_3.txt')
        elif name=='antipodal_neighbour':
            a=spins(d/'image_2.txt');b=spins(d/'image_3.txt');b[0][3:]=[-x for x in a[0][3:]];write_spins(d/'image_3.txt',b)
        elif name=='wrong_coordinates':
            rows=spins(d/'image_0.txt');rows[0][0]+=.2;write_spins(d/'image_0.txt',rows)
        elif name=='nonunit_spin':
            rows=spins(d/'image_0.txt');rows[0][3:]=[0.,0.,0.];write_spins(d/'image_0.txt',rows)
        elif name=='iteration_budget':(d/'input').write_text((d/'input').read_text().replace('gneb:max-iterations=200000','gneb:max-iterations=1'))
        r=subprocess.run([str(exe)],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60)
        assert r.returncode==expected,(name,r.returncode,r.stderr.decode(errors='replace'))
        if expected==2:assert not json.loads((d/'gneb_status.json').read_text())['numerically_converged']
        negative.append({'case':name,'expected_exit':expected,'actual_exit':r.returncode,'pass':True})
    # An existing output must be rejected, with its accepted status file unchanged.
    old=hashlib.sha256((out/'public9/gneb_status.json').read_bytes()).hexdigest()
    r=subprocess.run([str(exe)],cwd=out/'public9',stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60)
    assert r.returncode==1 and hashlib.sha256((out/'public9/gneb_status.json').read_bytes()).hexdigest()==old
    negative.append({'case':'output_overwrite','expected_exit':1,'actual_exit':1,'pass':True})
    # Coordinate-labelled rows may be permuted without changing the Hamiltonian.
    d=copy_case('permuted_static');rows=spins(d/'image_0.txt');write_spins(d/'permuted.txt',list(reversed(rows)))
    (d/'chain.txt').write_text('permuted.txt\n');(d/'input').write_text((d/'input').read_text().replace('gneb:mode=path','gneb:mode=evaluate'))
    r=subprocess.run([str(exe)],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60);assert r.returncode==0
    with (d/'gneb_energies.csv').open() as f:energy=float(next(csv.DictReader(f))['energy_meV'])
    assert abs(energy-result9['minimum_energy_meV'])<1e-8
    # Free endpoint minimization, then exact serial/parallel kernel comparison.
    relax_results=[]
    for nthreads in sorted(set([1,threads])):
        d=copy_case('relax_threads_'+str(nthreads));rows=spins(d/'image_0.txt')
        for i,row in enumerate(rows):
            angle=.015*math.sin(i*.3);x,y,z=row[3:]
            row[3:]=[math.cos(angle)*x+math.sin(angle)*z,y,-math.sin(angle)*x+math.cos(angle)*z]
        write_spins(d/'relax_seed.txt',rows);(d/'chain.txt').write_text('relax_seed.txt\n')
        (d/'input').write_text((d/'input').read_text().replace('gneb:mode=path','gneb:mode=relax').replace('gneb:force-tolerance=1e-6','gneb:force-tolerance=1e-7')+f'\ngneb:threads={nthreads}\ngneb:checkpoint-every=50\n')
        r=subprocess.run([str(exe)],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
        assert r.returncode==0,r.stderr.decode(errors='replace')
        status=json.loads((d/'gneb_status.json').read_text())
        with (d/'gneb_energies.csv').open() as f:en=float(next(csv.DictReader(f))['energy_meV'])
        assert status['numerically_converged'] and status['endpoint_torque_meV']<1e-7
        assert status['sampled_forward_barrier_meV'] is None and not status['path_barrier_evaluated']
        assert abs(en-result9['minimum_energy_meV'])<1e-8
        relax_results.append({'threads':nthreads,'energy_meV':en,'max_torque_meV':status['endpoint_torque_meV'],'iterations':status['iterations'],
                              'final_sha256':hashlib.sha256((d/'gneb_image_0.txt').read_bytes()).hexdigest()})
    assert len(set(x['final_sha256'] for x in relax_results))==1,'Parallel native relaxation differs from serial'
    if threads>1:
        d=copy_case('parallel_path');(d/'input').write_text((d/'input').read_text()+f'\ngneb:threads={threads}\n')
        r=subprocess.run([str(exe)],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
        assert r.returncode==0,r.stderr.decode(errors='replace')
        assert (d/'gneb_energies.csv').read_bytes()==(out/'public9/gneb_energies.csv').read_bytes()
        for k in range(9):assert (d/f'gneb_image_{k}.txt').read_bytes()==(out/f'public9/gneb_image_{k}.txt').read_bytes()
    d=copy_case('relax_budget');(d/'chain.txt').write_text('image_2.txt\n')
    (d/'input').write_text((d/'input').read_text().replace('gneb:mode=path','gneb:mode=relax').replace('gneb:max-iterations=200000','gneb:max-iterations=1'))
    r=subprocess.run([str(exe)],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60)
    assert r.returncode==2 and not json.loads((d/'gneb_status.json').read_text())['numerically_converged']
    negative.append({'case':'relax_iteration_budget','expected_exit':2,'actual_exit':2,'pass':True})
    result={'utc':datetime.now(timezone.utc).isoformat(),'all_passed':True,'data_scope':'Only public Spirit benchmark and explicitly defined analytic unit tests',
            'public_benchmarks':[result9,result17],'core_test_passed':True,'negative_tests':negative,'coordinate_permutation_passed':True,'endpoint_relaxation':relax_results,
            'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'platform':os.name,'compiler':subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]}
    report=Path(report);report.parent.mkdir(parents=True,exist_ok=True);report.write_text(json.dumps(result,indent=2)+'\n');print('All public acceptance tests passed.')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--exe',type=Path,required=True);p.add_argument('--cxx',default='g++');p.add_argument('--report',type=Path,default=ROOT/'build/acceptance.json');p.add_argument('--openmp-threads',type=int,default=1);a=p.parse_args();main(a.exe,a.cxx,a.report,a.openmp_threads)
