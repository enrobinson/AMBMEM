# AMBMEM

AMBMEM is a research-preview tool for constructing Amber-compatible membrane-only and protein-membrane systems.

## Validated capabilities

- POPC host membranes
- POPS and POPE reference support
- Cholesterol
- User-supplied single-residue custom sterols
- Membrane-only systems
- OPM-oriented protein-membrane systems
- OPC water
- KCl and NaCl
- Symmetric and asymmetric leaflet compositions
- Symmetric and asymmetric solvent buffers
- Reproducible seeded builds
- Exact per-leaflet composition
- Amber topology and restart output

## Requirements

AMBMEM must run inside a Conda environment containing:

- AmberTools 26 or newer
- Python 3.10 or newer
- NumPy
- PyYAML

## Installation from GitHub

```bash
conda create -n ambmem -c conda-forge \
    ambertools=26 python=3.12 numpy pyyaml pytest
conda activate ambmem

git clone https://github.com/enrobinson/AMBMEM.git
cd AMBMEM
python -m pip install .
```

For development:

```bash
python -m pip install -e .
```

Verify the environment:

```bash
ambmem doctor
```

## Quick start

```bash
cp examples/membrane_only_control.yaml build.yaml
ambmem validate build.yaml
ambmem build build.yaml
```

## Commands

```text
ambmem doctor
ambmem validate CONFIG.yaml
ambmem build CONFIG.yaml [--dry-run]
ambmem replicas CONFIG.yaml --seeds 1001 2003 3007 --output-prefix replica
ambmem validate-output topology.prmtop coordinates.rst7
ambmem inspect report.json
```

## Packaged references

AMBMEM packages validated coordinate references for POPC, POPS, POPE, and cholesterol. Amber force fields, Lipid21, OPC, ion parameters, GAFF2, ParmEd, and `tleap` come from the active AmberTools installation.

## Custom sterols

Custom parameters are not bundled. Users provide matching MOL2, Amber LIB/OFF, and FRCMOD files. See `docs/custom-sterols.md`.

## Documentation

- [Installation](docs/installation.md)
- [Quick start](docs/quickstart.md)
- [YAML reference](docs/yaml-reference.md)
- [Membrane-only builds](docs/membrane-only.md)
- [Protein-membrane builds](docs/protein-membrane.md)
- [Custom sterols](docs/custom-sterols.md)
- [Replicas](docs/replicas.md)
- [Outputs](docs/outputs.md)
- [Validation](docs/validation.md)
- [Troubleshooting](docs/troubleshooting.md)

## Status

Version 0.1.0 is an internal research preview. Additional Amber lipid types require characterization before full support is claimed.

## Simulation-readiness validation

AMBMEM validates the exact final `topology.prmtop` and `coordinates.rst7` after packing, sterol replacement, OPC placement, ion placement, topology generation, and periodic box assignment.

The final gate includes:

- exact topology/restart reload with ParmEd
- finite coordinate and box validation
- periodic minimum-image overlap scans
- hard failure for unrelated physical atoms below 0.5 A
- hard failure for nonbonded heavy atoms below 1.0 A
- CPPTRAJ overlap and abnormal-bond checks
- CPPTRAJ bonded and PME initial-energy evaluation
- hard failure for nonfinite energy or `|E| > 1e12 kJ/mol`

A passing output contains `final_validation.json`, CPPTRAJ diagnostics, and `SIMULATION_READY`. Do not minimize or simulate an output unless `SIMULATION_READY` is present.

