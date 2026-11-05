# Contributing

AMBMEM is currently an internal research preview.

1. Create a branch from `main`.
2. Add or update tests for every behavior change.
3. Run `python -m pytest`.
4. Run `ambmem doctor` and `ambmem validate` on an example configuration.
5. Do not commit protein coordinates, private custom parameters, trajectories, or generated Amber systems.
6. Document any newly validated lipid reference in the packaged reference-library manifest.
