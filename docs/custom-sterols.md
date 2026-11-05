# Custom sterols

Users provide a matching MOL2, Amber LIB/OFF, and FRCMOD. No project-specific custom parameters are bundled.

Requirements:

1. One complete residue in the MOL2 and library.
2. MOL2, library, YAML, and final residue names agree.
3. Atom names are unique and agree across files.
4. The residue charge matches `expected_formal_charge`; neutral sterols should sum to zero before replication.
5. `anchor_atom` identifies the interfacial polar atom.
6. `axis_end_atoms` define the distal membrane-core endpoint, optionally by centroid.
7. Placement uses proper rigid rotations; reflections are disabled by default.
8. Run `ambmem validate` to perform the MOL2 and tleap checks before packing.
