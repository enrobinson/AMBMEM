# Outputs

AMBMEM writes the completed simulation system, validation evidence, build
metadata, and several intermediate files to the configured output directory.

Not every file is produced in every build mode. Protein-membrane builds,
membrane-only builds, disabled ion placement, and failed builds may produce
slightly different subsets.

Most of these preliminary files do not need to be retained after a sucessful
packing, and are only needed to diagnose packing failures. 

## Primary simulation files

- `topology.prmtop`  
  Final Amber topology for the complete system. The topology includes the
  protein or other supplied solute, retained structural waters, membrane
  lipids and sterols as requested, generated solvent, and
  ions.

- `coordinates.rst7`  
  Final Amber restart coordinates corresponding to
  `topology.prmtop`. The restart also contains the final periodic box
  dimensions.

- `system_membrane_solvated.pdb`  
  PDB representation of the final system for visualization and inspection.
  The Amber `topology.prmtop` and `coordinates.rst7` pair remains the
  authoritative simulation input.

- `SIMULATION_READY`  
  Marker created only after the exact final topology-coordinate pair passes
  AMBMEM's final validation gate. The absence of this file means that the
  build did not complete validation successfully, and logs should be inspected.

## Build report

- `report.json`  
  
`report.json` is written only when the complete build reaches the reporting
  stage. A failed build may still contain useful earlier diagnostic files even
  when `report.json` is absent.

## Structural-water preprocessing

Protein-membrane builds that use an authoritative solute topology and
coordinates may write the following records.

- `structural_water_cleanup.json`  
  Audit of input structural waters. The report records which waters were
  accepted, retained with a warning, or rejected; the atom pairs and distances
  used for classification; the cleanup policy; the prepared solute files; and
  whether any rejected structural waters remained after preprocessing.

- `structural_water_blocks.json`  
  Describes the final atom blocks used to protect retained structural waters.
  The file records the fixed nonwater solute block, retained structural-water
  block, generated-environment block, atom and residue ranges, and whether
  topology-coordinate ordering was preserved.

- `authoritative_solute_preprocessing/`  
  Directory containing the topology-preserving preprocessing stages applied
  to the authoritative solute. Files in this directory provide an audit trail
  from the supplied PRMTOP/RST7 pair to the cleaned and reordered pair used in
  final assembly.

  The final prepared files normally include:

  - `03_clean_ordered.prmtop`  
    Prepared authoritative solute topology after structural-water cleanup and
    atom-order normalization.

  - `03_clean_ordered.rst7`  
    Coordinates corresponding exactly to `03_clean_ordered.prmtop`.

  Additional numbered files and CPPTRAJ scripts or logs may be present to
  document earlier preprocessing stages.

## Pre-ion environment files

- `environment_preion_lipids.pdb`  
  Lipid21 modular-lipid interchange PDB used to parameterize the membrane
  before ions are introduced. Standard phospholipids are written as complete
  modular lipids, such as `PA-PC-OL`, `PA-PE-OL`, or `PA-PS-OL`.

- `environment_preion_bulk.pdb`  
  Pre-ion PDB containing generated bulk molecules that are not written through
  the modular-lipid interchange file, principally solvent and other supported
  nonmodular species.

- `environment_preion.prmtop`  
  Amber topology of the generated membrane and solvent environment before ion
  replacement and before combination with the authoritative protein or solute
  topology.

- `environment_preion.rst7`  
  Coordinates and periodic box corresponding to
  `environment_preion.prmtop`.

- `tleap_environment_preion.log`  
  TLEaP output from parameterizing the pre-ion membrane and solvent
  environment. Consult this file for residue-template, atom-type, bonding, or
  parameterization warnings.

## Final environment files

- `environment_with_ions_lipids.pdb`  
  Lipid21 modular-lipid interchange PDB for the environment after selected
  solvent molecules have been replaced with ions.

- `environment_with_ions_bulk.pdb`  
  Bulk-molecule PDB for the final generated environment, including solvent and
  ions.

- `environment_with_ions.prmtop`  
  Amber topology of the generated environment after ion placement. In
  protein-membrane mode, this environment topology is subsequently combined
  with the authoritative solute topology.

- `environment_with_ions.rst7`  
  Coordinates corresponding to `environment_with_ions.prmtop`.

- `tleap_environment_final.log`  
  TLEaP output from parameterizing the final membrane, solvent, and ion
  environment.

## Final assembly intermediates

- `final_unboxed_intermediate.prmtop`  
  Intermediate combined topology created before the exact final periodic
  topology-coordinate pair is written.

- `final_unboxed_intermediate.rst7`  
  Coordinates corresponding to `final_unboxed_intermediate.prmtop`.

These files are retained to diagnose topology-combination, atom-ordering, or
periodic-box problems. They are not the recommended simulation inputs.

## Final validation

- `final_validation.json`  
  The report includes topology and coordinate paths, atom and residue counts, periodic box
  dimensions, validation thresholds, minimum-distance scans, CPPTRAJ results,
  initial energy, accepted pre-existing protein-contact warnings, fatal
  failures, and the final `simulation_ready` state.

- `final_validation.cpptraj`  
  CPPTRAJ input script generated by AMBMEM to validate the exact final
  `topology.prmtop` and `coordinates.rst7` pair.

- `final_validation_cpptraj.log`  
  Complete CPPTRAJ output from final structure and energy validation.

- `final_check_all_atoms.dat`  
  CPPTRAJ report of unrelated physical atom pairs closer than the final
  all-atom cutoff. A successful build normally produces an empty file.

- `final_check_heavy_atoms.dat`  
  CPPTRAJ report of nonbonded heavy-atom pairs closer than the final
  heavy-atom cutoff. A successful build normally produces an empty file.
  Protein-protein contacts already present in the authoritative input may be
  retained as explicit warnings when AMBMEM verifies that the same contacts
  existed before membrane construction.

- `final_check_bonds.dat`  
  CPPTRAJ structure report for abnormal bond lengths, very close atom pairs,
  and other supported structure checks. A successful build normally produces
  an empty file.

- `final_energy.dat`  
  CPPTRAJ initial energy decomposition for the exact final system. The output
  includes bond, angle, dihedral, 1-4, van der Waals, electrostatic, and total
  energy terms.

## Failure behavior

AMBMEM retains diagnostic files when a build fails whenever those files have
already been produced. For example:

- A failure during membrane placement may occur before topology files exist.
- A TLEaP failure may leave interchange PDB files and a TLEaP log.
- A final-validation failure normally leaves the exact final PRMTOP/RST7 pair,
  the CPPTRAJ script and log, the individual validation reports, and
  `final_validation.json`.
- `report.json` and `SIMULATION_READY` may be absent when the build terminates
  before successful completion.

Do not treat the existence of `topology.prmtop` and `coordinates.rst7` alone as
proof that a system is ready for simulation. A release-ready build must also
contain `SIMULATION_READY`, and `final_validation.json` must report:

```json
"simulation_ready": true
