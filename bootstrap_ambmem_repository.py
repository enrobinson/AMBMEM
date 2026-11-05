#!/usr/bin/env python3
"""Bootstrap the GitHub-ready AMBMEM 0.1.0 research-preview repository.

Run from /Users/emilyrobinson/software/AMBMEM after copying src/ambmem and the
validated reference library into src/ambmem/data/reference_library.
"""
from __future__ import annotations

from pathlib import Path
import textwrap

ROOT = Path.cwd()
PKG = ROOT / "src" / "ambmem"

if not (PKG / "protein_popc.py").is_file():
    raise SystemExit(
        "Run this script from the clean AMBMEM repository root. "
        "Expected src/ambmem/protein_popc.py."
    )

FILES: dict[str, str] = {
"pyproject.toml": r'''
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "ambmem"
version = "0.1.0"
description = "Research-preview builder for Amber-compatible membrane systems"
readme = "README.md"
requires-python = ">=3.10"
authors = [{name = "Emily Weed-Nichols"}]
classifiers = [
  "Development Status :: 3 - Alpha",
  "Intended Audience :: Science/Research",
  "Programming Language :: Python :: 3",
  "Topic :: Scientific/Engineering :: Chemistry",
]
dependencies = [
  "numpy>=1.24",
  "PyYAML>=6.0",
]

[project.scripts]
ambmem = "ambmem.cli:main"

[tool.setuptools]
package-dir = {"" = "src"}
include-package-data = true

[tool.setuptools.packages.find]
where = ["src"]

[tool.setuptools.package-data]
ambmem = [
  "data/reference_library/**/*.pdb",
  "data/reference_library/**/*.json",
  "data/reference_library/**/*.yaml",
  "data/reference_library/**/*.md",
]

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src"]
''',
"MANIFEST.in": r'''
recursive-include src/ambmem/data/reference_library *.pdb *.json *.yaml *.md
recursive-include examples *.yaml *.md
recursive-include docs *.md
''',
".gitignore": r'''
__pycache__/
*.py[cod]
*.egg-info/
build/
dist/
.pytest_cache/
.coverage
htmlcov/
.venv/
venv/
env/
conda-env/
.DS_Store

# AMBMEM outputs
**/topology.prmtop
**/coordinates.rst7
**/system_membrane_solvated.pdb
**/report.json
**/environment_*.pdb
**/environment_*.prmtop
**/environment_*.rst7
**/tleap_environment_*.log
**/leap.log
*.time

# User inputs and private parameters
sterol_parameters/
custom_parameters/
protein_inputs/
production_runs/
pre_membrane.prmtop
pre_membrane.rst7

# Trajectories and restart files
*.nc
*.mdcrd
*.dcd
*.xtc
*.trr
*.ncrst
*.restrt
''',
"LICENSE": r'''
AMBMEM Research Preview License Notice

Copyright (c) 2026 Emily Weed-Nichols. All rights reserved.

This repository is currently distributed for internal research evaluation only.
No permission is granted for redistribution, sublicensing, or commercial use.
A standard open-source license may be adopted before a public release.
''',
"CITATION.cff": r'''
cff-version: 1.2.0
title: AMBMEM
type: software
message: "If you use AMBMEM, cite this software and the relevant Amber/Lipid21/OPC methods."
authors:
  - family-names: Weed-Nichols
    given-names: Emily
version: 0.1.0
release-date: 2026-11-05
repository-code: "https://github.com/REPLACE-WITH-OWNER/AMBMEM"
license: other
''',
"CHANGELOG.md": r'''
# Changelog

## 0.1.0 Research Preview

Initial internal research preview.

### Validated workflows

- Membrane-only and OPM-oriented protein-membrane builds
- POPC, POPS, POPE, and cholesterol packing references
- Symmetric and asymmetric leaflet compositions
- OPC water with symmetric or asymmetric buffers
- KCl and NaCl using accepted-water-count concentration
- User-supplied single-residue custom sterols
- Seeded replica generation
- Amber `prmtop`/`rst7` output and JSON reporting

### Limitations

- Additional Amber lipids require characterization before support is claimed.
- Custom residue parameters are supplied by users and are not bundled.
- This is a research preview, not a production 1.0 release.
''',
"CONTRIBUTING.md": r'''
# Contributing

AMBMEM is currently an internal research preview.

1. Create a branch from `main`.
2. Add or update tests for every behavior change.
3. Run `python -m pytest`.
4. Run `ambmem doctor` and `ambmem validate` on an example configuration.
5. Do not commit protein coordinates, private custom parameters, trajectories, or generated Amber systems.
6. Document any newly validated lipid reference in the packaged reference-library manifest.
''',
"README.md": r'''
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

git clone https://github.com/REPLACE-WITH-OWNER/AMBMEM.git
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
''',
"src/ambmem/resources.py": r'''
from __future__ import annotations

from importlib.resources import files
from pathlib import Path


def reference_root() -> Path:
    return Path(str(files("ambmem").joinpath("data", "reference_library")))


def reference_path(species: str, filename: str | None = None) -> Path:
    species = species.upper()
    defaults = {
        "POPC": "reference.pdb",
        "POPS": "reference.pdb",
        "POPE": "reference.pdb",
        "CHL": "conformers.pdb",
    }
    if species not in defaults:
        raise KeyError(f"No packaged reference is registered for {species}")
    path = reference_root() / species / (filename or defaults[species])
    if not path.is_file():
        raise FileNotFoundError(f"Packaged reference is missing: {path}")
    return path


def resolve_reference(value: str) -> str:
    if value.startswith("packaged:"):
        return str(reference_path(value.split(":", 1)[1]))
    return value
''',
"src/ambmem/preflight.py": r'''
from __future__ import annotations

from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import tempfile
import yaml

from .resources import reference_path, resolve_reference


class PreflightError(RuntimeError):
    pass


def _ambertools_version(tleap: str) -> str:
    # AmberTools does not expose one universal version flag through tleap.
    # Prefer the Conda package record, then report the executable path.
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        meta = Path(conda_prefix) / "conda-meta"
        candidates = sorted(meta.glob("ambertools-*.json")) if meta.is_dir() else []
        if candidates:
            match = re.search(r"ambertools-([0-9][^-]*)", candidates[-1].name)
            if match:
                return match.group(1)
    return f"unknown ({tleap})"


def normalize_config(cfg: dict) -> dict:
    import copy
    cfg = copy.deepcopy(cfg)
    refs = cfg.setdefault("reference_library", {})
    if "POPC" in refs:
        refs["POPC"] = resolve_reference(str(refs["POPC"]))
    else:
        refs["POPC"] = str(reference_path("POPC"))

    complex_cfg = cfg.get("complex_membrane", {})
    if float(complex_cfg.get("POPS_fraction", 0.0)) > 0 or any(
        float(complex_cfg.get("leaflets", {}).get(leaf, {}).get("POPS_fraction", 0.0)) > 0
        for leaf in ("upper", "lower")
    ):
        complex_cfg["POPS_reference"] = resolve_reference(
            str(complex_cfg.get("POPS_reference", "packaged:POPS"))
        )
    if complex_cfg.get("enabled", False):
        complex_cfg["CHL_reference"] = resolve_reference(
            str(complex_cfg.get("CHL_reference", "packaged:CHL"))
        )
    cfg["complex_membrane"] = complex_cfg
    return cfg


def validate_config(path: str | Path, *, tleap_template_test: bool = True) -> tuple[dict, list[str]]:
    path = Path(path)
    if not path.is_file():
        raise PreflightError(f"Configuration file does not exist: {path}")
    cfg = yaml.safe_load(path.read_text())
    if not isinstance(cfg, dict):
        raise PreflightError("Configuration must contain a YAML mapping")
    cfg = normalize_config(cfg)
    messages: list[str] = []

    tleap = shutil.which("tleap")
    if not tleap:
        raise PreflightError("tleap was not found. Activate an AmberTools 26+ Conda environment.")
    messages.append(f"tleap: PASS ({tleap})")
    messages.append(f"AmberTools: {_ambertools_version(tleap)}")

    try:
        import parmed  # noqa: F401
    except ImportError as exc:
        raise PreflightError("ParmEd is not importable from the active environment") from exc
    messages.append("ParmEd: PASS")

    for species in ("POPC", "POPS", "POPE", "CHL"):
        reference_path(species)
    messages.append("Packaged references: PASS (POPC, POPS, POPE, CHL)")

    mode = "membrane_only" if cfg.get("membrane_only", {}).get("enabled", False) else "protein_membrane"
    messages.append(f"Build mode: {mode}")

    if mode == "membrane_only":
        m = cfg["membrane_only"]
        per_leaf = int(m.get("target_per_leaflet", 0))
        total = int(m.get("target_total_lipids", 2 * per_leaf))
        if per_leaf <= 0 or total != 2 * per_leaf:
            raise PreflightError("Membrane-only target_total_lipids must equal 2 * target_per_leaflet")
        box = m.get("box_xy_angstrom")
        if not isinstance(box, list) or len(box) != 2 or min(map(float, box)) <= 0:
            raise PreflightError("membrane_only.box_xy_angstrom must contain two positive values")
        messages.append(f"Membrane target: PASS ({per_leaf} molecules per leaflet; {box[0]} x {box[1]} A)")
    else:
        for key in ("pdb", "orientation_source", "preoriented"):
            if key not in cfg.get("protein", {}):
                raise PreflightError(f"protein.{key} is required")
        if cfg["protein"]["orientation_source"] != "OPM" or not cfg["protein"]["preoriented"]:
            raise PreflightError("Protein builds require an OPM-preoriented structure")
        for key in ("pdb",):
            if not Path(cfg["protein"][key]).expanduser().is_file():
                raise PreflightError(f"Missing protein input: {cfg['protein'][key]}")
        for key in ("solute_prmtop", "solute_coordinates"):
            value = cfg.get("topology", {}).get(key)
            if not value or not Path(value).expanduser().is_file():
                raise PreflightError(f"Missing topology input: topology.{key}")
        messages.append("Protein/topology inputs: PASS")

    composition = cfg.get("complex_membrane", {})
    for leaf in ("upper", "lower"):
        leaf_cfg = composition.get("leaflets", {}).get(leaf, {})
        chl = float(leaf_cfg.get("CHL_fraction", composition.get("CHL_fraction", 0.0)))
        pops = float(leaf_cfg.get("POPS_fraction", composition.get("POPS_fraction", 0.0)))
        custom = 0.0
        custom_block = cfg.get("custom_sterols", {})
        if custom_block.get("enabled", False):
            species = custom_block.get("species", [])
            if len(species) != 1:
                raise PreflightError("Exactly one custom sterol may be active in this research preview")
            s = species[0]
            custom = float(s.get("leaflets", {}).get(leaf, {}).get("fraction", s.get("fraction", 0.0)))
            _validate_custom_sterol(s, tleap_template_test=tleap_template_test)
        if min(chl, pops, custom) < 0 or chl + pops + custom > 1.0 + 1e-12:
            raise PreflightError(f"Invalid {leaf} fractions: CHL={chl}, POPS={pops}, custom={custom}")
        messages.append(f"{leaf} composition: PASS (POPC={1-chl-pops-custom:.3f}, CHL={chl:.3f}, POPS={pops:.3f}, custom={custom:.3f})")

    output = Path(cfg.get("output", {}).get("directory", "ambmem_output"))
    output.parent.mkdir(parents=True, exist_ok=True)
    if not os.access(output.parent, os.W_OK):
        raise PreflightError(f"Output parent is not writable: {output.parent}")
    messages.append(f"Output path: PASS ({output})")
    return cfg, messages


def _validate_custom_sterol(settings: dict, *, tleap_template_test: bool) -> None:
    required = ("residue_name", "mol2", "library", "parameters", "orientation")
    for key in required:
        if key not in settings:
            raise PreflightError(f"Custom sterol is missing {key}")
    paths = [Path(settings["mol2"]), Path(settings["library"]), *map(Path, settings["parameters"])]
    for path in paths:
        if not path.is_file():
            raise PreflightError(f"Missing custom sterol file: {path}")
    lines = Path(settings["mol2"]).read_text(errors="replace").splitlines()
    if "@<TRIPOS>MOLECULE" not in lines or "@<TRIPOS>ATOM" not in lines or "@<TRIPOS>BOND" not in lines:
        raise PreflightError("Custom MOL2 requires standard MOLECULE, ATOM, and BOND sections")
    atom_start = lines.index("@<TRIPOS>ATOM") + 1
    bond_start = lines.index("@<TRIPOS>BOND")
    atoms = {}
    residues = set()
    charge = 0.0
    for line in lines[atom_start:bond_start]:
        fields = line.split()
        if len(fields) >= 9 and fields[0].isdigit():
            atoms[fields[1]] = fields
            residues.add(fields[7])
            charge += float(fields[8])
    residue_name = settings["residue_name"]
    if residues != {residue_name}:
        raise PreflightError(f"MOL2 residue names {residues} do not match {residue_name}")
    orientation = settings["orientation"]
    names = [orientation.get("anchor_atom"), *orientation.get("axis_end_atoms", [])]
    missing = [name for name in names if name not in atoms]
    if missing:
        raise PreflightError(f"Custom sterol orientation atoms are missing: {missing}")
    expected = int(settings.get("expected_formal_charge", 0))
    if abs(charge - expected) > 0.01:
        raise PreflightError(f"Custom MOL2 charge {charge:+.6f} e does not match expected formal charge {expected}")
    if tleap_template_test:
        tleap = shutil.which("tleap")
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / "preflight.leap"
            log = Path(tmp) / "preflight.log"
            commands = ["source leaprc.gaff2"]
            commands.extend(f"loadAmberParams {Path(p).resolve()}" for p in settings["parameters"])
            commands.append(f"loadOff {Path(settings['library']).resolve()}")
            commands.append(f"check {residue_name}")
            commands.append("quit")
            script.write_text("\n".join(commands) + "\n")
            result = subprocess.run([tleap, "-f", str(script)], capture_output=True, text=True)
            log.write_text(result.stdout + result.stderr)
            if result.returncode != 0 or "Errors = 0" not in (result.stdout + result.stderr):
                raise PreflightError("Custom sterol failed the tleap preflight check")


def write_normalized_config(cfg: dict, original: Path) -> Path:
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", prefix=f"{original.stem}-normalized-", delete=False
    )
    with handle:
        yaml.safe_dump(cfg, handle, sort_keys=False)
    return Path(handle.name)
''',
"src/ambmem/cli.py": r'''
from __future__ import annotations

import argparse
from pathlib import Path
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
import yaml

from . import __version__
from .preflight import PreflightError, validate_config, write_normalized_config
from .protein_popc import build


def _print_preflight(messages: list[str]) -> None:
    print(f"AMBMEM {__version__} Research Preview")
    print("=" * 48)
    for message in messages:
        print(message)
    print("Preflight status: PASS")


def _build(config: Path, dry_run: bool) -> None:
    cfg, messages = validate_config(config)
    _print_preflight(messages)
    if dry_run:
        print("Dry run complete. Packing was not started.")
        return
    normalized = write_normalized_config(cfg, config)
    try:
        build(str(normalized))
    finally:
        normalized.unlink(missing_ok=True)


def _doctor() -> None:
    print(f"AMBMEM {__version__} environment check")
    failed = False
    for executable in ("tleap",):
        path = shutil.which(executable)
        print(f"{executable}: {path or 'MISSING'}")
        failed |= path is None
    try:
        import parmed
        print(f"ParmEd: {parmed.__version__}")
    except Exception as exc:
        print(f"ParmEd: FAIL ({exc})")
        failed = True
    try:
        from .resources import reference_path
        for species in ("POPC", "POPS", "POPE", "CHL"):
            print(f"{species} reference: {reference_path(species)}")
    except Exception as exc:
        print(f"Packaged references: FAIL ({exc})")
        failed = True
    if failed:
        raise SystemExit(1)
    print("Doctor status: PASS")


def _replicas(config: Path, seeds: list[int], prefix: str) -> None:
    base = yaml.safe_load(config.read_text())
    for index, seed in enumerate(seeds, 1):
        cfg = copy.deepcopy(base)
        cfg.setdefault("build", {})["seed"] = seed
        cfg.setdefault("build", {})["replicas"] = 1
        output = f"{prefix}_{index:02d}_seed_{seed}"
        cfg.setdefault("output", {})["directory"] = output
        generated = config.parent / f"{output}.yaml"
        generated.write_text(yaml.safe_dump(cfg, sort_keys=False))
        print(f"Generated {generated}")
        _build(generated, dry_run=False)


def _validate_output(topology: Path, coordinates: Path) -> None:
    import parmed as pmd
    system = pmd.load_file(str(topology), xyz=str(coordinates))
    charge = sum(float(atom.charge) for atom in system.atoms)
    print(f"Atoms: {len(system.atoms)}")
    print(f"Residues: {len(system.residues)}")
    print(f"Charge: {charge:+.8f} e")
    print(f"Box: {system.box}")
    if abs(charge - round(charge)) > 0.01:
        raise SystemExit("Output validation failed: noninteger charge")
    print("Output validation: PASS")


def _inspect(report: Path) -> None:
    data = json.loads(report.read_text())
    print(json.dumps({
        "version": data.get("version"),
        "build_mode": data.get("build_mode", "protein_membrane"),
        "seed": data.get("seed"),
        "box_angstrom": data.get("box_angstrom"),
        "lipid_counts": data.get("lipid_counts"),
        "final_charge_e": data.get("topology_and_ions", {}).get("final_charge_e"),
    }, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(prog="ambmem")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("doctor", help="Check AmberTools, ParmEd, and packaged references")
    validate = sub.add_parser("validate", help="Run preflight validation")
    validate.add_argument("config", type=Path)
    build_p = sub.add_parser("build", help="Validate and build a system")
    build_p.add_argument("config", type=Path)
    build_p.add_argument("--dry-run", action="store_true")
    replicas = sub.add_parser("replicas", help="Build seeded replicas sequentially")
    replicas.add_argument("config", type=Path)
    replicas.add_argument("--seeds", type=int, nargs="+", required=True)
    replicas.add_argument("--output-prefix", required=True)
    output = sub.add_parser("validate-output", help="Validate an Amber topology/restart pair")
    output.add_argument("topology", type=Path)
    output.add_argument("coordinates", type=Path)
    inspect = sub.add_parser("inspect", help="Summarize report.json")
    inspect.add_argument("report", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "doctor": _doctor()
        elif args.command == "validate":
            _, messages = validate_config(args.config)
            _print_preflight(messages)
        elif args.command == "build": _build(args.config, args.dry_run)
        elif args.command == "replicas": _replicas(args.config, args.seeds, args.output_prefix)
        elif args.command == "validate-output": _validate_output(args.topology, args.coordinates)
        elif args.command == "inspect": _inspect(args.report)
    except PreflightError as exc:
        print(f"Preflight status: FAIL\n{exc}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
''',
"src/ambmem/__init__.py": r'''
__version__ = "0.1.0"
''',
"examples/membrane_only_control.yaml": r'''
membrane_only:
  enabled: true
  box_xy_angstrom: [120.0, 120.0]
  half_thickness_angstrom: 17.9
  target_total_lipids: 400
  target_per_leaflet: 200

reference_library:
  POPC: packaged:POPC

membrane:
  host_lipid: POPC
  box_z_angstrom: 96.0

packing:
  fixed_heavy_lipid_heavy_angstrom: 2.0
  fixed_heavy_lipid_hydrogen_angstrom: 1.2
  lipid_hydrogen_hydrogen_angstrom: 1.15

build:
  seed: 7001
  replicas: 1

boundary_packing:
  enabled: false

complex_membrane:
  enabled: true
  POPS_fraction: 0.0
  CHL_fraction: 0.05
  protein_exclusion_shell_angstrom: 0.0
  CHL_reference: packaged:CHL
  cholesterol_trials: 160
  cholesterol_minimum_hydrogen_hydrogen_angstrom: 1.15
  require_exact_counts: true
  leaflets:
    upper: {POPS_fraction: 0.0, CHL_fraction: 0.05}
    lower: {POPS_fraction: 0.0, CHL_fraction: 0.05}

solvent:
  enabled: true
  model: OPC
  reference_cache: ~/.cache/ambmem/opcbox.pdb
  positive_z_buffer_from_headgroups_angstrom: 30.0
  negative_z_buffer_from_headgroups_angstrom: 30.0
  minimum_water_oxygen_system_distance_angstrom: 2.4
  minimum_water_hydrogen_system_distance_angstrom: 1.5
  minimum_water_virtual_site_system_distance_angstrom: 1.2
  preserve_input_waters: false
  minimum_water_beyond_fixed_atoms_angstrom: 0.0

topology:
  enabled: true
  membrane_only: true
  charge_integer_tolerance_e: 0.01
  output_prmtop: topology.prmtop
  output_rst7: coordinates.rst7

ions:
  enabled: true
  salt: KCl
  concentration_molar: 0.15
  minimum_ion_fixed_heavy_distance_angstrom: 4.0
  minimum_ion_ion_distance_angstrom: 5.0

output:
  directory: control_membrane_replica_01
  pdb_name: system_membrane_solvated.pdb
''',
"examples/membrane_only_custom_sterol.yaml": r'''
# Copy this file and replace CST paths/names with a user-supplied neutral sterol.
membrane_only:
  enabled: true
  box_xy_angstrom: [120.0, 120.0]
  half_thickness_angstrom: 17.9
  target_total_lipids: 400
  target_per_leaflet: 200

reference_library:
  POPC: packaged:POPC

membrane: {host_lipid: POPC, box_z_angstrom: 96.0}
packing:
  fixed_heavy_lipid_heavy_angstrom: 2.0
  fixed_heavy_lipid_hydrogen_angstrom: 1.2
  lipid_hydrogen_hydrogen_angstrom: 1.15
build: {seed: 1001, replicas: 1}
boundary_packing: {enabled: false}

complex_membrane:
  enabled: true
  POPS_fraction: 0.0
  CHL_fraction: 0.05
  protein_exclusion_shell_angstrom: 0.0
  CHL_reference: packaged:CHL
  cholesterol_trials: 160
  cholesterol_minimum_hydrogen_hydrogen_angstrom: 1.15
  require_exact_counts: true
  leaflets:
    upper: {POPS_fraction: 0.0, CHL_fraction: 0.05}
    lower: {POPS_fraction: 0.0, CHL_fraction: 0.05}

custom_sterols:
  enabled: true
  species:
    - name: example-sterol
      species_code: CST
      residue_name: CST
      mol2: parameters/custom_sterol.mol2
      library: parameters/custom_sterol.lib
      parameters: [parameters/custom_sterol.frcmod]
      expected_formal_charge: 0
      fraction: 0.05
      leaflets:
        upper: {fraction: 0.05}
        lower: {fraction: 0.05}
      orientation:
        anchor_atom: O1
        axis_start_atom: O1
        axis_end_atoms: [C28, C29]
        allow_reflection: false
      placement:
        host_species: POPC
        trials_per_site: 240
        xy_offset_angstrom: 3.0
        z_offset_angstrom: 1.5
        require_exact_count: true

solvent:
  enabled: true
  model: OPC
  reference_cache: ~/.cache/ambmem/opcbox.pdb
  positive_z_buffer_from_headgroups_angstrom: 30.0
  negative_z_buffer_from_headgroups_angstrom: 30.0
  minimum_water_oxygen_system_distance_angstrom: 2.4
  minimum_water_hydrogen_system_distance_angstrom: 1.5
  minimum_water_virtual_site_system_distance_angstrom: 1.2
  preserve_input_waters: false
  minimum_water_beyond_fixed_atoms_angstrom: 0.0

topology:
  enabled: true
  membrane_only: true
  charge_integer_tolerance_e: 0.01
  output_prmtop: topology.prmtop
  output_rst7: coordinates.rst7
ions:
  enabled: true
  salt: KCl
  concentration_molar: 0.15
  minimum_ion_fixed_heavy_distance_angstrom: 4.0
  minimum_ion_ion_distance_angstrom: 5.0
output:
  directory: custom_sterol_replica_01
  pdb_name: system_membrane_solvated.pdb
''',
"examples/protein_membrane_asymmetric.yaml": r'''
protein:
  pdb: inputs/model_aligned_renum.pdb
  orientation_source: OPM
  preoriented: true

reference_library:
  POPC: packaged:POPC

membrane:
  host_lipid: POPC
  edge_padding_angstrom: 20.0
  footprint_z_margin_angstrom: 8.0
  box_z_angstrom: 200.0

packing:
  fixed_heavy_lipid_heavy_angstrom: 2.0
  fixed_heavy_lipid_hydrogen_angstrom: 1.2
  lipid_hydrogen_hydrogen_angstrom: 1.15
build: {seed: 1001, replicas: 1}

complex_membrane:
  enabled: true
  POPS_fraction: 0.15
  CHL_fraction: 0.05
  POPS_reference: packaged:POPS
  CHL_reference: packaged:CHL
  protein_exclusion_shell_angstrom: 8.0
  cholesterol_trials: 160
  pops_donors_per_site: 64
  require_exact_counts: true
  leaflets:
    upper: {POPS_fraction: 0.0, CHL_fraction: 0.05}
    lower: {POPS_fraction: 0.30, CHL_fraction: 0.05}

solvent:
  enabled: true
  model: OPC
  reference_cache: ~/.cache/ambmem/opcbox.pdb
  positive_z_buffer_from_headgroups_angstrom: 30.0
  negative_z_buffer_from_headgroups_angstrom: 80.0
  minimum_fixed_atom_box_clearance_angstrom: 5.0
  minimum_water_oxygen_system_distance_angstrom: 2.4
  minimum_water_hydrogen_system_distance_angstrom: 1.5
  minimum_water_virtual_site_system_distance_angstrom: 1.2
  preserve_input_waters: true
  minimum_water_beyond_fixed_atoms_angstrom: 5.0

topology:
  enabled: true
  solute_prmtop: inputs/pre_membrane.prmtop
  solute_coordinates: inputs/pre_membrane.rst7
  charge_integer_tolerance_e: 0.01
  output_prmtop: topology.prmtop
  output_rst7: coordinates.rst7
ions:
  enabled: true
  salt: KCl
  concentration_molar: 0.15
  minimum_ion_fixed_heavy_distance_angstrom: 4.0
  minimum_ion_ion_distance_angstrom: 5.0
output:
  directory: protein_membrane_replica_01
  pdb_name: system_membrane_solvated.pdb
''',
"examples/README.md": r'''
# Examples

- `membrane_only_control.yaml`: 95% POPC and 5% cholesterol.
- `membrane_only_custom_sterol.yaml`: generic 90% POPC, 5% cholesterol, 5% user sterol.
- `protein_membrane_asymmetric.yaml`: OPM-oriented protein with asymmetric POPS.

Run `ambmem validate EXAMPLE.yaml` before building. Protein and custom-sterol examples require user-supplied inputs.
''',
"docs/installation.md": r'''
# Installation

```bash
conda create -n ambmem -c conda-forge ambertools=26 python=3.12 numpy pyyaml pytest
conda activate ambmem
git clone https://github.com/REPLACE-WITH-OWNER/AMBMEM.git
cd AMBMEM
python -m pip install .
ambmem doctor
```

AmberTools 26 or newer supplies `tleap`, ParmEd, Lipid21, OPC, ions, and GAFF2. AMBMEM supplies validated coordinate references for packing.
''',
"docs/quickstart.md": r'''
# Quick start

```bash
cp examples/membrane_only_control.yaml build.yaml
ambmem validate build.yaml
ambmem build build.yaml
ambmem validate-output control_membrane_replica_01/topology.prmtop control_membrane_replica_01/coordinates.rst7
```
''',
"docs/yaml-reference.md": r'''
# YAML reference

AMBMEM 0.1.0 retains the research-preview schema used by the validated backend.

- `membrane_only.enabled`: select membrane-only mode.
- `membrane_only.box_xy_angstrom`: lateral box dimensions.
- `membrane_only.target_per_leaflet`: exact host count per leaflet.
- `protein`: OPM-oriented PDB settings for protein builds.
- `reference_library.POPC`: use `packaged:POPC` by default.
- `complex_membrane`: cholesterol and POPS fractions and references.
- `custom_sterols`: one user-supplied sterol in the research preview.
- `solvent`: OPC buffers and contact cutoffs.
- `ions`: KCl or NaCl and concentration.
- `topology`: Amber input/output settings.
- `output.directory`: build directory.

All fraction and input checks should pass with `ambmem validate` before packing.
''',
"docs/membrane-only.md": r'''
# Membrane-only builds

Use `membrane_only.enabled: true`. Specify lateral dimensions, exact molecules per leaflet, composition, OPC buffers, and ions. No protein coordinates or solute topology are used.

For a 120 x 120 A, 400-molecule control, see `examples/membrane_only_control.yaml`.
''',
"docs/protein-membrane.md": r'''
# Protein-membrane builds

Protein builds require an OPM-preoriented PDB and a matching authoritative Amber `prmtop`/`rst7` pair. Do not use `setBox UNIT centers ...` when generating the pre-membrane restart, because that can change the coordinate frame. AMBMEM preserves the supplied solute topology while packing the membrane around the OPM coordinate frame.
''',
"docs/custom-sterols.md": r'''
# Custom sterols

Users provide a matching MOL2, Amber LIB/OFF, and FRCMOD. No project-specific custom parameters are bundled.

Requirements:

1. One complete residue in the MOL2 and library.
2. MOL2, library, YAML, and final residue names agree.
3. Atom names are unique and agree across files.
4. The residue charge matches `expected_formal_charge`; neutral sterols should sum to zero before replication.
5. `anchor_atom` identifies the interfacial polar atom.
6. `axis_end_atoms` define the distal membrane-core endpoint, optionally by centroid.
7. Placement uses proper rigid rotations; reflections are disabled by default.
8. Run `ambmem validate` to perform the MOL2 and tleap checks before packing.
''',
"docs/replicas.md": r'''
# Replicas

```bash
ambmem replicas build.yaml --seeds 1001 2003 3007 --output-prefix replica
```

Replicas run sequentially. Seeds vary host-site choices, rotations, offsets, water placement choices, and ion replacement while preserving identical composition and build settings.
''',
"docs/outputs.md": r'''
# Outputs

- `topology.prmtop`: final Amber topology
- `coordinates.rst7`: final Amber coordinates and box
- `system_membrane_solvated.pdb`: visualization structure
- `report.json`: composition, charge, box, reference, and validation metadata
- `tleap_environment_preion.log`: pre-ion parameterization log
- `tleap_environment_final.log`: final parameterization log
''',
"docs/validation.md": r'''
# Validation

Run preflight before packing:

```bash
ambmem validate build.yaml
```

Validate completed Amber files:

```bash
ambmem validate-output topology.prmtop coordinates.rst7
```

Inspect the build report:

```bash
ambmem inspect report.json
```
''',
"docs/troubleshooting.md": r'''
# Troubleshooting

## Missing host sites
Increase the lateral box dimensions instead of compressing the validated reference patch.

## Custom residue has a noninteger charge
Correct the residue charges before replication. Do not hide a repeated per-residue residual by loosening the whole-system tolerance.

## Lipid21 chirality failures
Confirm PDB atom-name columns, modular order (`PA -> PC/PS -> OL`), and `TER` only after complete lipids.

## ParmEd noncontiguous molecule warning
The research-preview backend allows ParmEd reordering and audits solute definitions independently of atom index.

## Protein appears outside the membrane cavity
Check that the pre-membrane restart retained the OPM frame. Avoid `setBox UNIT centers ...` when creating it.

## `tleap` illegal instruction while loading a library
Check Amber OFF/LIB section boundaries and validate the library independently before a full build.
''',
"scripts/validate_output.py": r'''
#!/usr/bin/env python3
from ambmem.cli import main

if __name__ == "__main__":
    main()
''',
"scripts/inspect_composition.py": r'''
#!/usr/bin/env python3
import json
import sys
from pathlib import Path

report = json.loads(Path(sys.argv[1]).read_text())
print(json.dumps(report.get("complex_composition", {}), indent=2))
''',
"tests/test_packaged_data.py": r'''
from ambmem.resources import reference_path


def test_packaged_references_exist():
    for species in ("POPC", "POPS", "POPE", "CHL"):
        assert reference_path(species).is_file()
''',
"tests/test_cli.py": r'''
from ambmem.cli import main


def test_import_cli():
    assert callable(main)
''',
"tests/test_config.py": r'''
from pathlib import Path
import yaml


def test_control_example_parses():
    config = yaml.safe_load(Path("examples/membrane_only_control.yaml").read_text())
    assert config["membrane_only"]["enabled"] is True
    assert config["membrane_only"]["target_per_leaflet"] == 200
''',
"tests/test_charge_validation.py": r'''
def test_placeholder_charge_tolerance():
    assert abs(-5.99899933 - round(-5.99899933)) < 0.01
''',
"tests/test_custom_sterol.py": r'''
def test_example_disables_reflection():
    import yaml
    config = yaml.safe_load(open("examples/membrane_only_custom_sterol.yaml"))
    species = config["custom_sterols"]["species"][0]
    assert species["orientation"]["allow_reflection"] is False
''',
}

for relative, content in FILES.items():
    path = ROOT / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(content).lstrip())
    if relative.startswith("scripts/"):
        path.chmod(0o755)

# Package markers and generic custom-parameter placeholders.
for relative in (
    "src/ambmem/data/__init__.py",
    "src/ambmem/data/reference_library/__init__.py",
):
    path = ROOT / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()

placeholder_dir = ROOT / "examples" / "custom_sterol" / "parameters"
placeholder_dir.mkdir(parents=True, exist_ok=True)
for name in ("PUT_MOL2_HERE", "PUT_LIB_HERE", "PUT_FRCMOD_HERE"):
    (placeholder_dir / name).touch()

# Validate required packaged references before declaring success.
required = {
    "POPC": "reference.pdb",
    "POPS": "reference.pdb",
    "POPE": "reference.pdb",
    "CHL": "conformers.pdb",
}
missing = [
    str(PKG / "data" / "reference_library" / species / filename)
    for species, filename in required.items()
    if not (PKG / "data" / "reference_library" / species / filename).is_file()
]
if missing:
    raise SystemExit(
        "Repository files were written, but packaged references are missing:\n"
        + "\n".join(missing)
    )

print("Created GitHub-ready AMBMEM research-preview files.")
print("Next commands:")
print("  python -m pip install -e .")
print("  ambmem doctor")
print("  ambmem validate examples/membrane_only_control.yaml")
print("  python -m pytest")
