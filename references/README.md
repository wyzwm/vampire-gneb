# Public reference inputs

The files under `spirit/` are unmodified source/test files from the public
[Spirit repository](https://github.com/spirit-code/spirit), commit
`e82250d3b14411c2c2fa292d143f13e3e111ad8c`.

|Local file|Original path|Purpose|
|---|---|---|
|`spirit/test_solvers.cpp`|`core/test/test_solvers.cpp`|Reference minimum and saddle energies|
|`spirit/solvers.cfg`|`core/test/input/solvers.cfg`|16×16 model parameters and boundaries|
|`spirit/fd_pairs.cfg`|`core/test/input/fd_pairs.cfg`|Explicit pair/D-vector convention|

Copyright and permission: [Spirit MIT license](../LICENSES/Spirit-MIT.txt).
Canonical LF SHA256 values are in [`upstream.json`](../upstream.json).
Generated benchmark spin images are in `examples/`; they are independently
reproduced from the public model and are not claimed to be upstream trajectory
files.
