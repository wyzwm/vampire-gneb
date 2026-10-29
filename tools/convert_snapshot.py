"""Convert an existing VAMPIRE coordinate/spin text pair; never launch a solver."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--coords',type=Path,required=True);p.add_argument('--spins',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
if a.out.exists() or Path(str(a.out)+'.json').exists():raise SystemExit('Output already exists')
coords=np.loadtxt(a.coords,skiprows=1)[:,2:5];s=np.loadtxt(a.spins,skiprows=1)
assert coords.shape==s.shape and s.ndim==2 and s.shape[1]==3 and np.isfinite(coords).all() and np.isfinite(s).all()
assert len(np.unique(np.round(coords,6),axis=0))==len(coords)
norm=np.linalg.norm(s,axis=1);error=float(np.max(abs(norm-1)))
assert error<=1e-5,'Spin norm differs from unit length beyond export tolerance'
# Keep original numbers. The native importer records and performs normalization.
with a.out.open('x',encoding='ascii') as f:
    f.write('GNEB_SPINS_V1\n%d\n'%len(s));np.savetxt(f,np.column_stack([coords,s]),fmt='%.17g')
record={'coords':str(a.coords.resolve()),'coords_SHA':hashlib.sha256(a.coords.read_bytes()).hexdigest(),
        'spins':str(a.spins.resolve()),'spins_SHA':hashlib.sha256(a.spins.read_bytes()).hexdigest(),
        'out_SHA':hashlib.sha256(a.out.read_bytes()).hexdigest(),'spin_count':len(s),'max_norm_error':error,
        'geometry_or_Hamiltonian_not_changed':True,'RNG_or_exact_continuation_not_claimed':True}
Path(str(a.out)+'.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2))
