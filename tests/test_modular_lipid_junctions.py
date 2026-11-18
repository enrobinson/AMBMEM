import numpy as np

from ambmem import protein_popc as builder


def _item(
    pa_c12,
    head_c11,
    head_c21,
    ol_c12,
    headgroup="PS",
    o12_distance=2.40,
):
    """Construct a minimal complete modular-lipid geometry fixture."""
    atoms=[
        {"resname":"PA","name":"C12","element":"C"},
        {"resname":headgroup,"name":"C11","element":"C"},
        {"resname":headgroup,"name":"O12","element":"O"},
        {"resname":headgroup,"name":"C21","element":"C"},
        {"resname":"OL","name":"C12","element":"C"},
    ]

    xyz=np.array([
        pa_c12,
        head_c11,
        [o12_distance,0.00,0.00],
        head_c21,
        ol_c12,
    ],dtype=float)

    return {
        "molecule":{"atoms":atoms},
        "xyz":xyz,
    }


def test_reference_like_pops_junctions_pass():
    item=_item(
        [0.00,0.00,0.00],
        [1.52,0.00,0.00],
        [4.00,0.00,0.00],
        [5.50,0.00,0.00],
    )

    valid,distances=builder._modular_lipid_junctions_valid(
        item,
        "PS",
        {},
    )

    assert valid
    assert 1.25 <= distances["PA_C12_to_headgroup_C11"] <= 1.75
    assert 1.25 <= distances["headgroup_C21_to_OL_C12"] <= 1.75


def test_stretched_pops_junction_is_rejected():
    item=_item(
        [0.00,0.00,0.00],
        [2.40,0.00,0.00],
        [4.00,0.00,0.00],
        [5.50,0.00,0.00],
    )

    valid,_=builder._modular_lipid_junctions_valid(
        item,
        "PS",
        {},
    )

    assert not valid


def test_compressed_pops_junction_is_rejected():
    item=_item(
        [0.00,0.00,0.00],
        [1.52,0.00,0.00],
        [4.00,0.00,0.00],
        [4.87,0.00,0.00],
    )

    valid,_=builder._modular_lipid_junctions_valid(
        item,
        "PS",
        {},
    )

    assert not valid


def _pops_item(o12_distance):
    atoms=[
        {"resname":"PA","name":"C12","element":"C"},
        {"resname":"PS","name":"C11","element":"C"},
        {"resname":"PS","name":"O12","element":"O"},
        {"resname":"PS","name":"C21","element":"C"},
        {"resname":"OL","name":"C12","element":"C"},
    ]

    return {
        "molecule":{"atoms":atoms},
        "xyz":np.array([
            [0.00,0.00,0.00],
            [1.52,0.00,0.00],
            [o12_distance,0.00,0.00],
            [4.00,0.00,0.00],
            [5.50,0.00,0.00],
        ],dtype=float),
    }


def test_pops_reference_like_internal_clearance_passes():
    item=_pops_item(2.40)

    valid,distances=(
        builder._modular_lipid_junctions_valid(
            item,
            "PS",
            {},
        )
    )

    assert valid
    assert distances["PA_C12_to_PS_O12"] == 2.40


def test_pops_observed_internal_clash_is_rejected():
    item=_pops_item(0.71)

    valid,distances=(
        builder._modular_lipid_junctions_valid(
            item,
            "PS",
            {},
        )
    )

    assert not valid
    assert distances["PA_C12_to_PS_O12"] < 1.75
