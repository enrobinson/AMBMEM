from ambmem.resources import reference_path


def test_packaged_references_exist():
    for species in ("POPC", "POPS", "POPE", "CHL"):
        assert reference_path(species).is_file()
