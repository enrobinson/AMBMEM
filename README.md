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
- ParmEd from the AmberTools environment
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

<!-- BEGIN AMBMEM-0.1.2-SAFETY -->
## Periodic membrane sizing

For membrane-only builds, `membrane_only.box_xy_angstrom` is an approximate lateral-size target. AMBMEM selects the nearest complete periodic supercell supported by the packaged POPC reference instead of cropping a noncommensurate box.

For example:

```yaml
membrane_only:
  enabled: true
  box_xy_angstrom: [120.0, 120.0]
```

The packaged POPC reference cell is approximately `63.491 x 63.805 A`. A 120 A target therefore selects a `2 x 2` supercell with an actual lateral box of approximately `126.982 x 127.610 A`.

AMBMEM retains all molecules in the complete supercell and converts requested leaflet percentages to the nearest realizable integer counts. For 256 molecules per leaflet:

- 95% POPC and 5% cholesterol becomes 243 POPC and 13 CHL per leaflet.
- 90% POPC, 5% cholesterol, and 5% custom sterol becomes 230 POPC, 13 CHL, and 13 custom sterol molecules per leaflet.

The requested percentages, selected repeats, actual box, integer counts, and realized fractions are recorded in `report.json`.

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

## Replica consistency

Replicas from the same composition and periodic supercell keep lipid, cholesterol, custom-sterol, and ion counts consistent and use deterministic dry atom ordering. Different seeds can change accepted OPC water counts because solvent-excluded volume changes slightly.

Users who require one shared topology may normalize waters afterward with ParmEd and validate that shared topology against every restart. AMBMEM does not apply hydrogen mass repartitioning.
<!-- END AMBMEM-0.1.2-SAFETY -->
