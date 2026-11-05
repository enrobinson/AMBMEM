from dataclasses import dataclass
from pathlib import Path
import numpy as np

@dataclass
class Atom:
    name: str
    resname: str
    resseq: int
    xyz: np.ndarray
    element: str

@dataclass
class Lipid:
    atoms: list
    anchor: int
    leaflet: str
    source_resseq: int


def read_pdb(path):
    """Read Lipid21-style whitespace PDB files; CRYST1 is optional."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f'Source PDB not found: {source}')

    atoms = []
    box = None

    for line_number, line in enumerate(
        source.read_text(errors='replace').splitlines(), start=1
    ):
        fields = line.split()

        if line.startswith('CRYST1'):
            if len(fields) < 4:
                raise ValueError(
                    f'Could not parse CRYST1 line {line_number}: {line}'
                )
            box = np.array(
                [float(fields[1]), float(fields[2]), float(fields[3])],
                dtype=float
            )
            continue

        if not line.startswith(('ATOM  ', 'HETATM')):
            continue

        if len(fields) < 8:
            raise ValueError(
                f'Could not parse atom line {line_number}: {line}'
            )

        try:
            name = fields[2]
            resname = fields[3]
            resseq = int(fields[4])
            xyz = np.array(
                [float(fields[5]), float(fields[6]), float(fields[7])],
                dtype=float
            )
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f'Could not parse atom line {line_number}: {line}'
            ) from exc

        element = next(
            (character.upper() for character in name if character.isalpha()),
            ''
        )

        atoms.append(
            Atom(
                name=name,
                resname=resname,
                resseq=resseq,
                xyz=xyz,
                element=element
            )
        )

    if not atoms:
        raise ValueError(f'No readable atoms found in {source}')

    return atoms, box


def group_lipids(atoms, component_resnames, anchor_name):
    """Group modular components with the same residue number into a lipid."""
    allowed = set(component_resnames)
    groups = {}

    for atom in atoms:
        if atom.resname in allowed:
            groups.setdefault(atom.resseq, []).append(atom)

    if not groups:
        raise ValueError(
            f'No components {sorted(allowed)} were found in the source file'
        )

    records = []
    anchor_z = []

    for resseq, group in groups.items():
        missing = allowed - {atom.resname for atom in group}
        if missing:
            raise ValueError(
                f'Molecule {resseq} is missing components {sorted(missing)}'
            )

        hits = [
            index for index, atom in enumerate(group)
            if atom.name == anchor_name
        ]
        if len(hits) != 1:
            raise ValueError(
                f'Molecule {resseq}: expected one {anchor_name}, '
                f'found {len(hits)}'
            )

        anchor = hits[0]
        anchor_z.append(group[anchor].xyz[2])
        records.append((resseq, group, anchor))

    sorted_z = np.sort(np.asarray(anchor_z, dtype=float))
    if len(sorted_z) == 1:
        midplane = 0.0
    else:
        midplane = float(
            (sorted_z[(len(sorted_z) - 1) // 2] +
             sorted_z[len(sorted_z) // 2]) / 2.0
        )

    lipids = []
    for resseq, group, anchor in records:
        leaflet = (
            'upper' if group[anchor].xyz[2] > midplane else 'lower'
        )
        lipids.append(
            Lipid(
                atoms=group,
                anchor=anchor,
                leaflet=leaflet,
                source_resseq=resseq
            )
        )

    return lipids, midplane
