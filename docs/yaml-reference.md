# YAML reference

AMBMEM 0.1.0 retains the research-preview schema used by the validated backend.

- `membrane_only.enabled`: select membrane-only mode.
- `membrane_only.box_xy_angstrom`: lateral box dimensions.
- `membrane_only.target_per_leaflet`: exact host count per leaflet.
- `protein`: OPM-oriented PDB settings for protein builds.
- `reference_library.POPC`: use `packaged:POPC` by default.
- `complex_membrane`: cholesterol and POPS fractions and references.
- `custom_sterols`: one user-supplied sterol in the research preview.
- `solvent`: OPC buffers and contact cutoffs.
- `ions`: KCl or NaCl and concentration.
- `topology`: Amber input/output settings.
- `output.directory`: build directory.

All fraction and input checks should pass with `ambmem validate` before packing.
