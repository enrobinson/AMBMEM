from pathlib import Path

import pytest
import ambmem.protein_popc as builder


def test_packaged_popc_cryst1_dimensions():
    path = Path(builder.__file__).parent / "data/reference_library/POPC/reference.pdb"
    assert path.is_file()
    cryst1 = next(line for line in path.read_text().splitlines() if line.startswith("CRYST1"))
    assert float(cryst1[6:15]) == pytest.approx(63.491, abs=0.001)
    assert float(cryst1[15:24]) == pytest.approx(63.805, abs=0.001)


def test_120_angstrom_target_selects_two_by_two_supercell():
    rx = max(1, round(120.0 / 63.491))
    ry = max(1, round(120.0 / 63.805))
    assert (rx, ry) == (2, 2)
    assert rx * 63.491 == pytest.approx(126.982, abs=0.001)
    assert ry * 63.805 == pytest.approx(127.610, abs=0.001)


def test_realisable_leaflet_compositions():
    total = 256
    chl = round(total * 0.05)
    custom = round(total * 0.05)
    assert (total - chl, chl) == (243, 13)
    assert (total - chl - custom, chl, custom) == (230, 13, 13)
