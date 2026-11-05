# Protein-membrane builds

Protein builds require an OPM-preoriented PDB and a matching authoritative Amber `prmtop`/`rst7` pair. Do not use `setBox UNIT centers ...` when generating the pre-membrane restart, because that can change the coordinate frame. AMBMEM preserves the supplied solute topology while packing the membrane around the OPM coordinate frame.
