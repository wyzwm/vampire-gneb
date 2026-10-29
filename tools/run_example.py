"""Run the shipped public Spirit benchmark in a new output directory."""
from pathlib import Path
import argparse,csv,hashlib,json,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1]
def verify(directory):
    directory=Path(directory);status=json.loads((directory/'gneb_status.json').read_text())
    reference=json.loads((directory/'provenance.json').read_text())
    with (directory/'gneb_energies.csv').open() as f:rows=list(csv.DictReader(f))
    energies=[float(r['energy_meV']) for r in rows]
    minimum=energies[0];saddle=max(energies);tolerance=reference['absolute_energy_tolerance_meV']
    if not status['numerically_converged']:raise AssertionError('GNEB did not converge')
    if abs(minimum-reference['reference_minimum_energy_meV'])>=tolerance:raise AssertionError('Minimum disagrees with public reference')
    if abs(saddle-reference['reference_saddle_energy_meV'])>=tolerance:raise AssertionError('Saddle disagrees with public reference')
    if abs(status['sampled_forward_barrier_meV']-(saddle-minimum))>1e-8:raise AssertionError('Barrier accounting mismatch')
    if status['max_GNEB_force_meV']>=1e-6:raise AssertionError('Force criterion not met')
    image=directory/f"gneb_image_{status['climbing_image']}.txt"
    spinlines=image.read_text().splitlines()[2:]
    moment_z=2*sum(float(line.split()[5]) for line in spinlines)/len(spinlines)
    result={'public_source_commit':reference['source_commit'],'images':len(rows),'spins_per_image':status['spins_per_image'],
            'minimum_energy_meV':minimum,'saddle_energy_meV':saddle,'forward_barrier_meV':saddle-minimum,
            'reference_minimum_energy_meV':reference['reference_minimum_energy_meV'],'reference_saddle_energy_meV':reference['reference_saddle_energy_meV'],
            'reference_forward_barrier_meV':reference['reference_forward_barrier_meV'],'absolute_energy_tolerance_meV':tolerance,
            'maximum_force_meV_per_radian':status['max_GNEB_force_meV'],'iterations':status['iterations'],
            'saddle_magnetization_z_muB_per_spin':moment_z,'reference_saddle_magnetization_z_muB_per_spin':2*.96657,
            'energy_and_force_checks_passed':True}
    (directory/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def run(executable,output):
    executable=Path(executable).resolve();output=Path(output).resolve();source=ROOT/'examples/spirit-skyrmion-16x16'
    if output.exists():raise SystemExit('Output directory exists; use a new directory')
    reference=json.loads((source/'provenance.json').read_text())
    for name,digest in reference['files'].items():
        if hashlib.sha256((source/name).read_bytes()).hexdigest()!=digest:raise SystemExit('Example checksum mismatch: '+name)
    output.mkdir(parents=True)
    for name in list(reference['files'])+['provenance.json']:shutil.copyfile(source/name,output/name)
    r=subprocess.run([str(executable)],cwd=output,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
    (output/'stdout.log').write_bytes(r.stdout);(output/'stderr.log').write_bytes(r.stderr)
    if r.returncode:raise SystemExit('Solver exited '+str(r.returncode)+'; see output logs')
    result=verify(output);print(json.dumps(result,indent=2));return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--exe',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.exe,a.output)
