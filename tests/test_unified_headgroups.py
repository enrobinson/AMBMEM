from ambmem.protein_popc import _allocate_leaflet_counts


def test_benchmark_outer_counts():
    assert _allocate_leaflet_counts(256,{"POPC":.60,"POPE":.10,"POPS":0,"CHL":.30}) == {"POPC":154,"POPE":25,"POPS":0,"CHL":77}

def test_benchmark_inner_counts():
    assert _allocate_leaflet_counts(256,{"POPC":.30,"POPE":.25,"POPS":.20,"CHL":.25}) == {"POPC":77,"POPE":64,"POPS":51,"CHL":64}

def test_counts_sum_exactly():
    result=_allocate_leaflet_counts(256,{"POPC":.60,"POPE":.10,"POPS":0,"CHL":.30})
    assert sum(result.values())==256
