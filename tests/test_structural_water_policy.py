import numpy as np
import parmed as pmd

from ambmem.protein_popc import _structural_water_clash_audit


def _system(distance):
    system=pmd.Structure()
    protein=pmd.Atom(name='C',type='C',atomic_number=6,mass=12.01,charge=0.0)
    system.add_atom(protein,'ALA',1)
    for name,number,mass in [('O',8,16.0),('H1',1,1.008),('H2',1,1.008),('EPW',0,0.0)]:
        atom=pmd.Atom(name=name,type=name,atomic_number=number,mass=mass,charge=0.0)
        system.add_atom(atom,'WAT',2)
    system.coordinates=np.array([
        [0.,0.,0.],[distance,0.,0.],[distance,.8,0.],
        [distance,-.8,0.],[distance,.1,0.],
    ])
    return system


def test_remove_policy_removes_complete_hard_clashing_water():
    cleaned,audit=_structural_water_clash_audit(
        _system(.5),
        {'structural_waters':{'hard_clash_policy':'remove'}},
    )
    assert audit['removed_structural_waters']==1
    assert audit['retained_structural_waters']==0
    assert len(cleaned.atoms)==1


def test_warning_water_is_retained():
    cleaned,audit=_structural_water_clash_audit(
        _system(1.2),
        {'structural_waters':{
            'hard_clash_policy':'remove',
            'heavy_atom_cutoff_angstrom':1.0,
            'warning_cutoff_angstrom':1.5,
        }},
    )
    assert audit['warning_retained_structural_waters']==1
    assert audit['removed_structural_waters']==0
    assert len(cleaned.atoms)==5


def test_accepted_water_is_retained_and_counted_for_ions():
    cleaned,audit=_structural_water_clash_audit(
        _system(2.0),
        {'structural_waters':{'hard_clash_policy':'remove'}},
    )
    assert audit['accepted_structural_waters']==1
    assert audit['counted_in_total_water_for_ion_target'] is True
    assert audit['ion_replacement_eligible'] is False
