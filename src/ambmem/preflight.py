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
