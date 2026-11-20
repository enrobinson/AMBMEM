# Final-system validation

AMBMEM 0.1.2 validates the exact final Amber topology/restart pair after every construction stage. The validator uses ParmEd, NumPy, and CPPTRAJ from AmberTools 26 or newer.

## Validation sequence

1. Reload `topology.prmtop` and `coordinates.rst7`.
2. Verify matching atom counts, finite coordinates, and positive box dimensions.
3. Scan periodic minimum-image distances while excluding direct bonds and virtual sites such as OPC `EPW`.
4. Run CPPTRAJ all-atom overlap, heavy-atom overlap, and abnormal-bond checks.
5. Evaluate bonded and PME energy for the initial frame with CPPTRAJ.
6. Write diagnostics and create `SIMULATION_READY` only if every check passes.

## Construction safeguards

AMBMEM builds complete commensurate supercells from the reference `CRYST1` vectors. Periodic minimum-image checks are applied during cholesterol and custom-sterol replacement, OPC placement, and ion placement. Every accepted OPC water is registered before later waters are tested.


