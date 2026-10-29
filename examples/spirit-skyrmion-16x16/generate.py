"""Reproduce the public Spirit 16x16 solver-test input; requires its helper binary.
All parameters come from references/spirit/solvers.cfg. No private material data.
"""
from pathlib import Path
import argparse, hashlib, json, subprocess

MEV_J=1.602176634e-22

def prepare(out, generator):
    out=Path(out).resolve();generator=Path(generator).resolve()
    if out.exists():raise SystemExit('Output directory exists; choose a new directory')
    out.mkdir(parents=True)
    r=subprocess.run([str(generator),'--prepare-only'],cwd=out,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=120)
    (out/'preparation_stdout.log').write_bytes(r.stdout);(out/'preparation_stderr.log').write_bytes(r.stderr)
    if r.returncode:raise SystemExit('Public benchmark endpoint preparation failed')
    lines=['1 1 1','1 0 0','0 1 0','0 0 1','1 1','0 .5 .5 .5 0 1 0','4 tensorial']
    for index,(dx,dy) in enumerate([(1,0),(0,1),(-1,0),(0,-1)]):
        j=10*MEV_J;Dx=6*dx*MEV_J;Dy=6*dy*MEV_J
        tensor=[j,0,-Dy,0,j,Dx,Dy,-Dx,j]
        lines.append(f'{index} 0 0 {dx} {dy} 0 '+' '.join(f'{v:.17g}' for v in tensor))
    (out/'file.ucf').write_text('\n'.join(lines)+'\n',encoding='ascii')
    (out/'Co.mat').write_text('''# Generic public square-lattice spin model, not a material parameterization.
material:num-materials=1
material[1]:material-name=PublicBenchmarkSpin
material[1]:damping-constant=1
material[1]:initial-spin-direction=0,0,1
material[1]:atomic-spin-moment=2!muB
material[1]:uniaxial-anisotropy-constant=0
material[1]:uniaxial-anisotropy-direction=0,0,1
''',encoding='ascii')
    (out/'input').write_text('''# Spirit official 16x16 skyrmion-collapse benchmark, serial GNEB.
create:crystal-structure=sc
create:periodic-boundaries-x
create:periodic-boundaries-y
create:spin-initialisation-random-seed=22
dimensions:system-size-x=1.6!nm
dimensions:system-size-y=1.6!nm
dimensions:system-size-z=0.1!nm
material:file=Co.mat
material:unit-cell-file=file.ucf
sim:program=gneb
sim:integrator=llg-heun
sim:temperature=0
sim:applied-field-strength=25!T
sim:applied-field-unit-vector=0,0,1
gneb:mode=path
gneb:chain-file=chain.txt
gneb:output-prefix=gneb
gneb:spring-constant=1
gneb:force-tolerance=1e-6
gneb:endpoint-tolerance=1e-6
gneb:max-iterations=200000
gneb:climb-after=2000
gneb:initial-step=0.005
gneb:maximum-step=0.05
gneb:report-every=500
''',encoding='ascii')
    # Canonical LF bytes keep the recorded hashes reproducible across Git checkouts.
    for file in out.iterdir():
        if file.is_file() and file.suffix!='.log':
            content=file.read_text(encoding='utf-8')
            with file.open('w',encoding='utf-8',newline='\n') as stream:stream.write(content)
    provenance={'source_repository':'https://github.com/spirit-code/spirit','source_commit':'e82250d3b14411c2c2fa292d143f13e3e111ad8c',
        'source_files':['core/test/test_solvers.cpp','core/test/input/solvers.cfg','core/src/utility/Configurations.cpp'],
        'source_parameters':{'grid':[16,16,1],'boundary_conditions':[1,1,0],'J_meV':10,'D_meV':6,'D_chirality':'Bloch +1','mu_muB':2,'Bz_T':25,'K_meV':0,'images':9,'initial_skyrmion_radius':5,'initial_phase_degrees':-90},
        'reference_minimum_energy_meV':-5849.69140625,'reference_saddle_energy_meV':-5811.5244140625,
        'reference_forward_barrier_meV':38.1669921875,'absolute_energy_tolerance_meV':3e-4,
        'tolerance_reason':'Published test values are float32; accounts for their print/representation rounding and documented Bohr-magneton constant difference.',
        'spin_images':'Independently generated from the public model and configuration recipe, not downloaded reference trajectories.',
        'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.suffix!='.log'}}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--generator',type=Path,required=True)
    args=parser.parse_args();prepare(args.out,args.generator)
