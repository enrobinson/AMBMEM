import copy
import numpy as np
import parmed as pmd

from ambmem.protein_popc import _split_fixed_solute_water_blocks


def _test_structure():
    structure=pmd.Structure()
    atoms=[]
    for name,resname,atomic_number,mass in [
        ('CA','ALA',6,12.01),
        ('O','WAT',8,16.00),
        ('H1','WAT',1,1.008),
        ('H2','WAT',1,1.008),
        ('EPW','WAT',0,0.0),
    ]:
        atom=pmd.Atom(name=name,type=name,atomic_number=atomic_number,mass=mass,charge=0.0)
        structure.add_atom(atom,resname,1 if resname=='ALA' else 2)
        atoms.append(atom)
    structure.coordinates=np.array([
        [0.0,0.0,0.0],[1.0,0.0,0.0],[1.8,0.0,0.0],
        [0.8,0.7,0.0],[1.0,0.2,0.0],
    ])
    return structure


def test_structural_water_split_preserves_counts_and_coordinates():
    source=_test_structure()
    nonwater,waters,audit=_split_fixed_solute_water_blocks(source)
    assert len(nonwater.atoms)==1
    assert len(waters.atoms)==4
    assert len(waters.residues)==1
    assert audit['retained_structural_water_residues']==1
    assert audit['retained_structural_waters_strippable'] is False
    assert np.allclose(waters.coordinates,source.coordinates[1:])


def test_source_structure_is_not_modified():
    source=_test_structure()
    original=copy.deepcopy(source.coordinates)
    _split_fixed_solute_water_blocks(source)
    assert len(source.atoms)==5
    assert np.allclose(source.coordinates,original)
