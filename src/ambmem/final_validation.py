from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import json
import math
import re
import shutil
import subprocess
import tempfile
import numpy as np


class FinalValidationError(RuntimeError):
    pass


def _is_virtual(atom) -> bool:
    return (
        float(getattr(atom, "mass", 0.0) or 0.0) < 0.1
        or int(getattr(atom, "atomic_number", 0) or 0) == 0
        or atom.name.strip().upper() in {"EP", "EPW", "LP", "LP1", "LP2", "MW", "VS"}
    )


def _label(atom) -> str:
    residue = atom.residue
    return f"{atom.idx + 1}:{residue.name.strip()}:{residue.idx + 1}:{atom.name.strip()}"


def _minimum_image(delta: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    return delta - lengths * np.rint(delta / lengths)


def _periodic_pair_scan(structure, unrelated_cutoff: float, heavy_cutoff: float, report_limit: int):
    xyz = np.asarray(structure.coordinates, dtype=float).reshape((-1, 3))
    box = np.asarray(structure.box, dtype=float)
    if box.size < 6 or np.any(~np.isfinite(box)) or np.any(box[:3] <= 0):
        raise FinalValidationError("Final restart has no valid periodic box")
    if np.max(np.abs(box[3:6] - 90.0)) > 1.0e-4:
        raise FinalValidationError("Research-preview validator currently requires an orthorhombic box")
    lengths = box[:3]
    physical = [i for i, atom in enumerate(structure.atoms) if not _is_virtual(atom)]
    bonded = {tuple(sorted((bond.atom1.idx, bond.atom2.idx))) for bond in structure.bonds}

    width = max(unrelated_cutoff, heavy_cutoff, 1.0)
    ncell = np.maximum(1, np.floor(lengths / width).astype(int))
    wrapped = xyz - lengths * np.floor((xyz + 0.5 * lengths) / lengths)
    cells = defaultdict(list)
    for index in physical:
        key = tuple(np.floor((wrapped[index] + 0.5 * lengths) / lengths * ncell).astype(int) % ncell)
        cells[key].append(index)

    severe, heavy = [], []
    minimum_all = (math.inf, None)
    minimum_heavy = (math.inf, None)
    visited = set()
    offsets = [(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)]
    for key, indices in cells.items():
        for offset in offsets:
            neighbor = tuple((np.asarray(key) + np.asarray(offset)) % ncell)
            for i in indices:
                for j in cells.get(neighbor, ()):
                    pair_index = (i, j) if i < j else (j, i)
                    if i == j or pair_index in visited:
                        continue
                    visited.add(pair_index)
                    if pair_index in bonded:
                        continue
                    delta = _minimum_image(xyz[j] - xyz[i], lengths)
                    distance = float(np.linalg.norm(delta))
                    pair = {
                        "distance_angstrom": distance,
                        "atom1": _label(structure.atoms[i]),
                        "atom2": _label(structure.atoms[j]),
                    }
                    if distance < minimum_all[0]:
                        minimum_all = (distance, pair)
                    if distance < unrelated_cutoff:
                        severe.append(pair)
                    atom_i, atom_j = structure.atoms[i], structure.atoms[j]
                    is_heavy = (
                        int(getattr(atom_i, "atomic_number", 0) or 0) > 1
                        and int(getattr(atom_j, "atomic_number", 0) or 0) > 1
                    )
                    if is_heavy:
                        if distance < minimum_heavy[0]:
                            minimum_heavy = (distance, pair)
                        if distance < heavy_cutoff:
                            heavy.append(pair)

    severe.sort(key=lambda item: item["distance_angstrom"])
    heavy.sort(key=lambda item: item["distance_angstrom"])
    return {
        "minimum_unrelated_distance_angstrom": None if minimum_all[1] is None else minimum_all[0],
        "minimum_unrelated_pair": minimum_all[1],
        "minimum_nonbonded_heavy_distance_angstrom": None if minimum_heavy[1] is None else minimum_heavy[0],
        "minimum_nonbonded_heavy_pair": minimum_heavy[1],
        "unrelated_pairs_below_cutoff": len(severe),
        "heavy_pairs_below_cutoff": len(heavy),
        "worst_unrelated_pairs": severe[:report_limit],
        "worst_heavy_pairs": heavy[:report_limit],
    }


def _parse_cpptraj_energy(path: Path):
    if not path.is_file():
        raise FinalValidationError("CPPTRAJ did not write final_energy.dat")
    headers = []
    values = None
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            headers = stripped.lstrip("#").split()
            continue
        parts = stripped.split()
        try:
            row = [float(value) for value in parts]
        except ValueError:
            continue
        if len(row) >= 2:
            values = row
            break
    if values is None:
        raise FinalValidationError("CPPTRAJ energy output has no numeric data row")
    # CPPTRAJ energy data are in kcal/mol. The first column is the frame index.
    component_values = values[1:]
    names = headers[-len(values):][1:] if len(headers) >= len(values) else [f"component_{i+1}" for i in range(len(component_values))]
    components_kcal = dict(zip(names, component_values))
    total_kcal = float(sum(component_values))
    return {
        "components_kcal_mol": components_kcal,
        "total_kcal_mol": total_kcal,
        "total_kj_mol": total_kcal * 4.184,
    }


def _run_cpptraj(prmtop: Path, rst7: Path, workdir: Path, bond_offset: float, cutoff: float):
    cpptraj = shutil.which("cpptraj")
    if not cpptraj:
        raise FinalValidationError("cpptraj was not found in the active AmberTools environment")
    input_path = workdir / "final_validation.cpptraj"
    log_path = workdir / "final_validation_cpptraj.log"
    all_report = workdir / "final_check_all_atoms.dat"
    heavy_report = workdir / "final_check_heavy_atoms.dat"
    bond_report = workdir / "final_check_bonds.dat"
    energy_path = workdir / "final_energy.dat"
    input_path.write_text(
        f"parm {prmtop.resolve()}\n"
        f"trajin {rst7.resolve()} 1 1\n"
        f"check !@EPW,EP,LP,LP1,LP2,MW,VS cut 0.5 nobondcheck reportfile {all_report.resolve()}\n"
        f"check !@H=,EPW,EP,LP,LP1,LP2,MW,VS cut 1.0 nobondcheck reportfile {heavy_report.resolve()}\n"
        f"check !@EPW,EP,LP,LP1,LP2,MW,VS cut 0.5 offset {bond_offset:.6f} reportfile {bond_report.resolve()}\n"
        f"energy InitialEnergy bond angle dihedral nb14 elec vdw etype pme cut {cutoff:.6f} out {energy_path.resolve()}\n"
        "run\n"
        "quit\n"
    )
    result = subprocess.run([cpptraj, "-i", str(input_path)], capture_output=True, text=True)
    log_text = result.stdout + "\n" + result.stderr
    log_path.write_text(log_text)
    if result.returncode != 0 or re.search(r"(^|\n)Error:", log_text, re.IGNORECASE):
        raise FinalValidationError(f"CPPTRAJ validation failed; see {log_path}")
    energy = _parse_cpptraj_energy(energy_path)

    def problem_lines(path):
        if not path.is_file():
            return []
        return [
            line.strip() for line in path.read_text(errors="replace").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]

    all_problems = problem_lines(all_report)
    heavy_problems = problem_lines(heavy_report)
    bond_problems = problem_lines(bond_report)
    return {
        "executable": cpptraj,
        "returncode": result.returncode,
        "log": str(log_path),
        "all_atom_report": str(all_report),
        "heavy_atom_report": str(heavy_report),
        "bond_report": str(bond_report),
        "energy_report": str(energy_path),
        "all_atom_problem_count": len(all_problems),
        "heavy_atom_problem_count": len(heavy_problems),
        "bond_problem_count": len(bond_problems),
        "all_atom_problems": all_problems[:100],
        "heavy_atom_problems": heavy_problems[:100],
        "bond_problems": bond_problems[:100],
        "energy": energy,
    }


def validate_final_system(
    prmtop: str | Path,
    rst7: str | Path,
    *,
    unrelated_cutoff_angstrom: float = 0.5,
    heavy_cutoff_angstrom: float = 1.0,
    bond_offset_angstrom: float = 0.5,
    energy_threshold_kj_mol: float = 1.0e12,
    cpptraj_cutoff_angstrom: float = 10.0,
    report_limit: int = 100,
    report_path: str | Path | None = None,
):
    import parmed as pmd

    prmtop, rst7 = Path(prmtop), Path(rst7)
    if not prmtop.is_file() or not rst7.is_file():
        raise FinalValidationError(f"Missing final Amber files: {prmtop}, {rst7}")
    structure = pmd.load_file(str(prmtop), xyz=str(rst7))
    xyz = np.asarray(structure.coordinates, dtype=float).reshape((-1, 3))
    if len(structure.atoms) != len(xyz):
        raise FinalValidationError("Final topology and restart atom counts do not match")
    if np.any(~np.isfinite(xyz)):
        raise FinalValidationError("Final coordinates contain NaN or Inf")

    pair_report = _periodic_pair_scan(structure, unrelated_cutoff_angstrom, heavy_cutoff_angstrom, report_limit)
    cpptraj_report = _run_cpptraj(prmtop, rst7, prmtop.parent, bond_offset_angstrom, cpptraj_cutoff_angstrom)
    energy_kj = cpptraj_report["energy"]["total_kj_mol"]
    failures = []
    if pair_report["unrelated_pairs_below_cutoff"]:
        failures.append(f"{pair_report['unrelated_pairs_below_cutoff']} unrelated atom pairs below {unrelated_cutoff_angstrom} A")
    if pair_report["heavy_pairs_below_cutoff"]:
        failures.append(f"{pair_report['heavy_pairs_below_cutoff']} nonbonded heavy pairs below {heavy_cutoff_angstrom} A")
    if cpptraj_report["all_atom_problem_count"]:
        failures.append(f"CPPTRAJ reported {cpptraj_report['all_atom_problem_count']} all-atom overlap/bond problems")
    if cpptraj_report["heavy_atom_problem_count"]:
        failures.append(f"CPPTRAJ reported {cpptraj_report['heavy_atom_problem_count']} heavy-atom overlap problems")
    if cpptraj_report["bond_problem_count"]:
        failures.append(f"CPPTRAJ reported {cpptraj_report['bond_problem_count']} overlap or abnormal-bond problems")
    if not np.isfinite(energy_kj) or abs(energy_kj) > energy_threshold_kj_mol:
        failures.append(f"CPPTRAJ initial energy {energy_kj:.8g} kJ/mol failed threshold {energy_threshold_kj_mol:.8g}")

    report = {
        "validator": "AMBMEM exact-final-pair CPPTRAJ gate",
        "prmtop": str(prmtop.resolve()),
        "rst7": str(rst7.resolve()),
        "atoms": len(structure.atoms),
        "residues": len(structure.residues),
        "box": None if structure.box is None else [float(value) for value in structure.box],
        "periodic_minimum_image": True,
        "thresholds": {
            "unrelated_atom_angstrom": unrelated_cutoff_angstrom,
            "nonbonded_heavy_atom_angstrom": heavy_cutoff_angstrom,
            "bond_offset_angstrom": bond_offset_angstrom,
            "absolute_energy_kj_mol": energy_threshold_kj_mol,
        },
        "pair_scan": pair_report,
        "cpptraj": cpptraj_report,
        "bond_ring_intersection": {
            "implemented": False,
            "status": "pending_for_public_1.0",
            "note": "Catastrophic distance, bond-length, and CPPTRAJ energy gates are active.",
        },
        "simulation_ready": not failures,
        "failures": failures,
    }
    report_path = Path(report_path) if report_path else prmtop.parent / "final_validation.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    marker = prmtop.parent / "SIMULATION_READY"
    marker.unlink(missing_ok=True)
    if failures:
        raise FinalValidationError("; ".join(failures))
    marker.write_text("AMBMEM exact final prmtop/rst7 validation passed.\n")
    return report
