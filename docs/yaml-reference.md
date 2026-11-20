# YAML reference


- `membrane_only.enabled`: select membrane-only mode.
- `membrane_only.box_xy_angstrom`: lateral box dimensions.
- `protein`: OPM-oriented PDB settings for protein builds.
- `custom_sterols`: one user-supplied sterol in the research preview.
- `solvent`: OPC buffers.
- `ions`: KCl or NaCl and concentration.
- `topology`: Amber input/output settings.
- `output.directory`: build directory.

All fraction and input checks should pass with `ambmem validate` before packing.
For full YAML example file see examples/
