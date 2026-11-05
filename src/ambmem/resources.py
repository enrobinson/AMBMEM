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
