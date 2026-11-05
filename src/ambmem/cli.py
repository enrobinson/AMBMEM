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
