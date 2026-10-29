# Free endpoints and optional CPU threading

Version 0.2 adds native single-image FIRE minimization, using the same VAMPIRE
energy and gradient as the path solver. Supply exactly one coordinate-labelled
image in `chain.txt`, with the normal geometry and material input:

```text
sim:program=gneb
sim:temperature=0
gneb:mode=relax
gneb:chain-file=chain.txt
gneb:force-tolerance=1e-7
gneb:max-iterations=200000
gneb:report-every=50
gneb:checkpoint-every=1000
```

Use a fresh directory for every endpoint. The final state is
`gneb_image_0.txt`; `gneb_status.json` reports the residual torque and convergence.
No physical time is advanced. No spin is pinned or frozen. An iteration budget
exit is code 2 and is **not** convergence. A converged gradient is necessary,
but does not certify a local minimum: inspect topology, finite perturbations and
curvature as appropriate before using a state as a GNEB endpoint. A disappearing
texture during free relaxation must not be reported as a metastable endpoint.

For periodic systems, use the complete periodic geometry at both endpoints and
every path image. Isolating a defect by cutting out a patch changes the boundary
value problem. Clearly identify which defects disappear and which survive.

## Optional OpenMP native kernel

```bash
python3 tools/build.py --source build/vampire --jobs 2 --openmp
python3 tests/run_tests.py --exe build/bin/vampire-gneb --openmp-threads 2
```

Set `gneb:threads=2` (or the number of exclusively allocated physical CPUs).
This runs **one process**, not MPI. VAMPIRE native field and site-energy loops
are divided into nonoverlapping atom ranges. Global energy is summed in the
same atom order as the single-thread run. The path optimizer and image loop
remain sequential, so scaling is not proportional to the thread count.

Parallel support is deliberately restricted to bilinear exchange (including
tensor DMI), second-order uniaxial anisotropy and static Zeeman energy. Other
anisotropies are rejected when `threads > 1`. The default is one thread; a
non-OpenMP build rejects a larger value. Do not use `mpirun` to launch duplicate
copies. The scheduler/cpuset must assign the requested physical CPUs, and the
operator must configure thread affinity, e.g. `OMP_PROC_BIND=true` and an
appropriate `GOMP_CPU_AFFINITY` for GNU OpenMP. Confirm actual thread bindings.

The public acceptance suite compares the native 16x16 endpoint relaxation and
nine-image GNEB result between one and multiple threads, including final image
bytes. It also tests the relaxed endpoint energy against the public Spirit
reference. This is implementation validation, not a universal scaling claim.

FIRE: Bitzek et al., Phys. Rev. Lett. 97, 170201 (2006),
https://doi.org/10.1103/PhysRevLett.97.170201.
Spin-manifold path method: Bessarab et al., CPC 196, 335 (2015),
https://doi.org/10.1016/j.cpc.2015.07.001.
