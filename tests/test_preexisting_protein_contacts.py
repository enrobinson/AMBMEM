from ambmem import final_validation as validator


class Residue:
    def __init__(self, name, idx):
        self.name = name
        self.idx = idx


class Atom:
    def __init__(self, residue, name):
        self.residue = residue
        self.name = name


def test_stable_pair_identity_is_order_independent():
    glu = Residue("GLU", 236)
    asn = Residue("ASN", 700)

    atom1 = Atom(glu, "OE2")
    atom2 = Atom(asn, "OD1")

    forward = validator._stable_pair_identity(atom1, atom2)
    reverse = validator._stable_pair_identity(atom2, atom1)

    assert forward == reverse
    assert forward == (
        ("ASN", 701, "OD1"),
        ("GLU", 237, "OE2"),
    )


def test_cpptraj_pair_identity_matches_stable_identity():
    line = (
        "1 Warning: Atoms 3748:GLU_237@OE2 and "
        "10873:ASN_701@OD1 are close (0.74)"
    )

    observed = validator._cpptraj_pair_identity(line)

    assert observed == (
        ("ASN", 701, "OD1"),
        ("GLU", 237, "OE2"),
    )


def test_non_contact_cpptraj_line_has_no_pair_identity():
    assert (
        validator._cpptraj_pair_identity(
            "1 Warning: Unusual bond length"
        )
        is None
    )


def test_protein_residue_classification():
    assert validator._is_protein_atom(
        Atom(Residue("GLU", 0), "OE2")
    )
    assert not validator._is_protein_atom(
        Atom(Residue("POPC", 0), "C1")
    )
