from pathlib import Path

import ambmem.protein_popc as builder


def test_unified_converter_constructs_headgroup_between_tails():
    source = Path(builder.__file__).read_text()

    sequence = """atoms=(
        retained_by_module['PA']['atoms']
        + transformed_atoms
        + retained_by_module['OL']['atoms']
    )"""

    assert sequence in source


def test_writer_audit_is_headgroup_aware():
    source = Path(builder.__file__).read_text()

    assert "expected=['PA',headgroup,'OL']" in source


def test_module_order_report_includes_pope():
    source = Path(builder.__file__).read_text()

    assert "'POPE':['PA','PE','OL']" in source
