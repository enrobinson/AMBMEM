from copy import deepcopy

from ambmem.protein_popc import (
    _canonicalize_environment_item_order,
)


def item(species, leaflet, site):
    return {
        "species": species,
        "leaflet": leaflet,
        "canonical_site_id": site,
    }


def signatures(items):
    return [
        (
            entry["species"],
            entry["leaflet"],
            entry["canonical_site_id"],
        )
        for entry in items
    ]


def test_canonical_order_is_seed_independent():
    replica_a = [
        item("FFM", "lower", (1, 1, 4)),
        item("POPC", "upper", (0, 0, 2)),
        item("CHL", "upper", (1, 0, 1)),
        item("FFM", "upper", (0, 1, 3)),
        item("POPC", "lower", (0, 0, 5)),
    ]

    replica_b = [
        item("POPC", "lower", (0, 0, 5)),
        item("FFM", "upper", (0, 1, 3)),
        item("FFM", "lower", (1, 1, 4)),
        item("CHL", "upper", (1, 0, 1)),
        item("POPC", "upper", (0, 0, 2)),
    ]

    ordered_a, report_a = (
        _canonicalize_environment_item_order(
            deepcopy(replica_a)
        )
    )

    ordered_b, report_b = (
        _canonicalize_environment_item_order(
            deepcopy(replica_b)
        )
    )

    assert signatures(ordered_a) == signatures(ordered_b)
    assert report_a["seed_independent_topology_order"]
    assert report_b["seed_independent_topology_order"]


def test_cholesterol_precedes_custom_sterol():
    items = [
        item("TMS", "upper", (0, 0, 2)),
        item("CHL", "upper", (0, 0, 1)),
        item("POPC", "upper", (0, 0, 3)),
    ]

    ordered, _ = _canonicalize_environment_item_order(
        items
    )

    assert [
        entry["species"]
        for entry in ordered
    ] == [
        "POPC",
        "CHL",
        "TMS",
    ]
