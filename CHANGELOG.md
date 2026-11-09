# Changelog

## 0.1.3 Research Preview

### Canonical replica topology ordering

- Assigns every periodic host-lipid site a stable canonical site ID.
- Preserves the site ID when POPC is replaced by cholesterol, POPS,
  POPE, or a user-supplied custom sterol.
- Sorts complete membrane molecules deterministically before solvation
  and Amber topology generation.
- Uses seed-independent ordering by species, leaflet, and host-site ID.
- Maintains identical dry residue and atom ordering across replicas
  generated from the same composition.
- Maintains identical ion counts and ion ordering across replicas.
- Allows OPC water counts to vary with solvent-exclusion geometry.
- Supports one shared topology after complete-water normalization.

### Internal validation

- Three FFM replicas had identical dry and ion ordering.
- After matching all replicas to 14,211 complete OPC waters, one
  canonical FFM topology passed the full periodic and CPPTRAJ gate
  with all three restart files.
- TMS replicas were tested using the same canonical-ordering and
  matched-water workflow.
- AMBMEM does not perform hydrogen mass repartitioning. HMR may be
  applied once afterward to the validated shared topology.

## 0.1.2 Research Preview

### Critical fixes

- Build complete commensurate periodic membrane supercells from packaged reference `CRYST1` dimensions.
- Treat requested membrane X/Y dimensions as approximate targets.
- Apply periodic minimum-image checks during sterol, OPC, and ion placement.
- Register accepted OPC waters immediately in the collision structure.
- Validate the exact written Amber topology/restart pair with periodic scans and CPPTRAJ overlap, bond, and energy checks.
- Create `SIMULATION_READY` only after all final checks pass.

### Internal validation

Three independent replicas passed the final gate for each of three conditions: 95% POPC plus 5% cholesterol; 90% POPC plus 5% cholesterol plus 5% user-supplied custom sterol; and a second 90/5/5 custom-sterol condition. Internal custom parameters are not distributed with AMBMEM.

## 0.1.0 Research Preview

Initial internal research preview.

### Validated workflows

- Membrane-only and OPM-oriented protein-membrane builds
- POPC, POPS, POPE, and cholesterol packing references
- Symmetric and asymmetric leaflet compositions
- OPC water with symmetric or asymmetric buffers
- KCl and NaCl using accepted-water-count concentration
- User-supplied single-residue custom sterols
- Seeded replica generation
- Amber `prmtop`/`rst7` output and JSON reporting

### Limitations

- Additional Amber lipids require characterization before support is claimed.
- Custom residue parameters are supplied by users and are not bundled.
- This is a research preview, not a production 1.0 release.
