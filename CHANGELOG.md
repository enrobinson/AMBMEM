# Changelog

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
