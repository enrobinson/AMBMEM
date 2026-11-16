import numpy as np

from ambmem.protein_popc import (
    _cholesterol_spread_candidate_order,
    _allocate_leaflet_counts,
)


def _grid_items(n=16):
    side=int(np.sqrt(n))
    return [
        {
            "anchor": np.array([float(i % side), float(i // side), 0.0]),
            "species": "POPC",
            "leaflet": "upper",
            "canonical_site_id": (0, 0, i),
        }
        for i in range(n)
    ]


def test_spread_order_is_complete_and_reproducible():
    items=_grid_items()
    cfg={"membrane_only":{"box_xy_angstrom":[4.0,4.0]}}
    a=_cholesterol_spread_candidate_order(list(range(16)),items,cfg,np.random.default_rng(7001))
    b=_cholesterol_spread_candidate_order(list(range(16)),items,cfg,np.random.default_rng(7001))
    assert a==b
    assert sorted(a)==list(range(16))


def test_high_cholesterol_benchmark_target_is_exact():
    result=_allocate_leaflet_counts(
        256,
        {"POPC":0.60,"POPE":0.10,"POPS":0.0,"CHL":0.30},
    )
    assert result=={"POPC":154,"POPE":25,"POPS":0,"CHL":77}
    assert sum(result.values())==256
