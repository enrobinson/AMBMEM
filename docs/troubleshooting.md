# Troubleshooting

## Missing host sites
Increase the lateral box dimensions instead of compressing the validated reference patch.

## Custom residue has a noninteger charge
Correct the residue charges before replication. Do not hide a repeated per-residue residual by loosening the whole-system tolerance.

## Lipid21 chirality failures
Confirm PDB atom-name columns, modular order (`PA -> PC/PS -> OL`), and `TER` only after complete lipids.

## ParmEd noncontiguous molecule warning
The research-preview backend allows ParmEd reordering and audits solute definitions independently of atom index.

## Protein appears outside the membrane cavity
Check that the pre-membrane restart retained the OPM frame. Avoid `setBox UNIT centers ...` when creating it.

## `tleap` illegal instruction while loading a library
Check Amber OFF/LIB section boundaries and validate the library independently before a full build.
