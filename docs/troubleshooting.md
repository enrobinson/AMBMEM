# Troubleshooting

## No SIMULATION_READY marker

The output failed final validation and must not be simulated. Inspect `final_validation.json`, `final_validation_cpptraj.log`, `final_check_all_atoms.dat`, `final_check_heavy_atoms.dat`, `final_check_bonds.dat`, and `final_energy.dat`.

## Actual box differs from requested box

This is expected. AMBMEM maps the requested X/Y target to the nearest complete periodic supercell defined by the packaged reference `CRYST1` dimensions.

## Replica topology hashes differ

Compare dry atom order and molecule counts. OPC water counts can differ across seeds. If dry atom order is identical, waters may be normalized later with ParmEd for a shared-topology workflow.

## Custom residue is nonneutral

Correct the MOL2 and Amber LIB charges so the residue charge matches `expected_formal_charge`. A small per-residue residual accumulates when many copies are inserted.

## Lipid21 connection errors

AMBMEM writes Lipid21 modules in validated order and places `TER` only after each complete lipid. Do not place `TER` between PA, PC/PS/PE, and OL modules.
