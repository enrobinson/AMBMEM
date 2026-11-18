import inspect

import ambmem.protein_popc as builder


def test_writer_recognizes_pope_as_modular_lipid():
    source = inspect.getsource(
        builder._write_lipid21_environment_files
    )

    compact = source.replace(" ", "")

    assert "{'POPC','POPE','POPS'}" in compact
    assert "'POPE':'PE'" in compact


def test_writer_uses_species_headgroup_mapping():
    source = inspect.getsource(
        builder._write_lipid21_environment_files
    )

    assert "headgroup_by_species" in source
    assert "headgroup:1" in source.replace(" ", "")


def test_converter_generates_correct_module_sequence():
    source = inspect.getsource(
        builder._convert_popc_headgroup_item
    )

    pa = source.index(
        "retained_by_module['PA']['atoms']"
    )

    headgroup = source.index(
        "+ transformed_atoms"
    )

    ol = source.index(
        "+ retained_by_module['OL']['atoms']"
    )

    assert pa < headgroup < ol
