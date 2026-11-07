# Final-system validation

AMBMEM 0.1.2 validates the exact final Amber topology/restart pair after every construction stage. The validator uses ParmEd, NumPy, and CPPTRAJ from AmberTools 26 or newer. OpenMM and sander are not required.

## Default configuration

```yaml
final_validation:
  unrelated_atom_cutoff_angstrom: 0.5
  nonbonded_heavy_cutoff_angstrom: 1.0
  bond_length_offset_angstrom: 0.5
  cpptraj_cutoff_angstrom: 10.0
  energy_threshold_kj_mol: 1.0e12
  report_limit: 100
```

These are catastrophic-assembly thresholds, not equilibration criteria.

## Validation sequence

1. Reload `topology.prmtop` and `coordinates.rst7`.
2. Verify matching atom counts, finite coordinates, and positive box dimensions.
3. Scan periodic minimum-image distances while excluding direct bonds and virtual sites such as OPC `EPW`.
4. Run CPPTRAJ all-atom overlap, heavy-atom overlap, and abnormal-bond checks.
5. Evaluate bonded and PME energy for the initial frame with CPPTRAJ.
6. Write diagnostics and create `SIMULATION_READY` only if every check passes.

## Construction safeguards

AMBMEM builds complete commensurate supercells from the reference `CRYST1` vectors. Periodic minimum-image checks are applied during cholesterol and custom-sterol replacement, OPC placement, and ion placement. Every accepted OPC water is registered before later waters are tested.

## Replica water counts

Seeded replicas can accept different OPC water counts because sterol coordinates change solvent exclusion. Dry composition and dry atom order remain consistent. Users who require one topology can normalize waters with ParmEd and then validate the shared topology against every restart.
