# Attribution, licenses and scope

This repository is an **unofficial extension** for VAMPIRE 7.0.0. It is not
maintained or endorsed by the VAMPIRE or Spirit authors.

## VAMPIRE

- Upstream: <https://github.com/richard-evans/vampire>
- Base tag: `7.0.0`
- Commit: `525bc27ee44c525aee229570f30f3d4c61d54f66`
- Preserve the upstream `license`, `BSD_licence`, per-file notices, and bundled
  third-party notices. The GPL version 2 text is retained as this repository's
  `LICENSE`; upstream BSD text is in `LICENSES/VAMPIRE-BSD.txt`.

New extension code is provided under GPL-2.0-or-later. This statement does not
relicense any upstream or third-party file. The extension reuses VAMPIRE's
initialized energy and field routines and must be distributed together with
their applicable notices and corresponding source when distributing a build.

## Spirit

- Upstream: <https://github.com/spirit-code/spirit>
- Reference commit: `e82250d3b14411c2c2fa292d143f13e3e111ad8c`
- Reference implementation: `core/src/engine/Method_GNEB.cpp` and
  `core/src/engine/Manifoldmath.cpp`.
- Public benchmark: `core/test/test_solvers.cpp` and
  `core/test/input/solvers.cfg`.
- Copyright (c) 2015 Gideon Müller and collaborators. MIT license, retained in
  full at `LICENSES/Spirit-MIT.txt`.

The energy-weighted tangent, product-sphere projection, spring and climbing-image
force definitions informed this implementation. A dependency-free kernel and
FIRE optimizer with sphere exponential updates and parallel transport were
written for this extension. No Spirit binary is linked. The benchmark files in
`references/spirit/` are unmodified MIT-licensed upstream files; generated spin
images are separately identified as reproductions from that public model.

## Integration changes beyond the new program

The patch also adds the missing standard `<cstdint>` include, corrects Qhull's
pointer-width type on Win64/MinGW, and bounds the anisotropy parameter buffer by
the actual material count before accessing the magnetic-moment array. These
changes fix compilation and array-bounds issues; they do not introduce new
interaction terms.

## Methods to cite

- Evans et al., *J. Phys.: Condens. Matter* **26**, 103202 (2014),
  <https://doi.org/10.1088/0953-8984/26/10/103202> — VAMPIRE/spin modeling.
- Müller et al., *Physical Review B* **99**, 224414 (2019),
  <https://doi.org/10.1103/PhysRevB.99.224414> — Spirit framework.
- Bessarab, Uzdin and Jónsson, *Computer Physics Communications* **196**, 335–347
  (2015), <https://doi.org/10.1016/j.cpc.2015.07.001> — magnetic GNEB.
- Bitzek et al., *Physical Review Letters* **97**, 170201 (2006),
  <https://doi.org/10.1103/PhysRevLett.97.170201> — FIRE.
- Henkelman and Jónsson, *J. Chem. Phys.* **113**, 9978 (2000),
  <https://doi.org/10.1063/1.1323224> — improved tangent.
- Henkelman, Uberuaga and Jónsson, *J. Chem. Phys.* **113**, 9901 (2000),
  <https://doi.org/10.1063/1.1329672> — climbing image.

Passing the software benchmarks is not a certification of a material's energy
barrier, stability, or lifetime. Version 0.1 is serial CPU only.
