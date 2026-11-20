"""Internal AMBMEM protocol constants.

These values define the validated AMBMEM construction and safety
protocol. Users should modify these values with caution, as
changes to the validated safety checks may lead to clashes in the
final produced structure.
"""
# -------------------------------------------------------------------
# Exact-final-pair validation
# -------------------------------------------------------------------
FINAL_UNRELATED_ATOM_CUTOFF_ANGSTROM = 0.5
FINAL_NONBONDED_HEAVY_CUTOFF_ANGSTROM = 1.0
FINAL_BOND_LENGTH_OFFSET_ANGSTROM = 0.5
FINAL_CPPTRAJ_CUTOFF_ANGSTROM = 10.0
FINAL_ENERGY_THRESHOLD_KJ_MOL = 1.0e12
FINAL_VALIDATION_REPORT_LIMIT = 100
# -------------------------------------------------------------------
# Charge and topology
# -------------------------------------------------------------------
CHARGE_INTEGER_TOLERANCE_E = 0.01
# -------------------------------------------------------------------
# Structural and generated water handling
# -------------------------------------------------------------------
MINIMUM_WATER_OXYGEN_SYSTEM_DISTANCE_ANGSTROM = 2.4
MINIMUM_WATER_HYDROGEN_SYSTEM_DISTANCE_ANGSTROM = 1.5
MINIMUM_WATER_VIRTUAL_SITE_SYSTEM_DISTANCE_ANGSTROM = 1.2
MINIMUM_WATER_BEYOND_FIXED_ATOMS_ANGSTROM = 0.0
PRESERVE_INPUT_WATERS = True
OPC_REFERENCE_CACHE = "~/.cache/ambmem/opcbox.pdb"
# -------------------------------------------------------------------
# Packaged references and host membrane
# -------------------------------------------------------------------
HOST_LIPID = "POPC"
POPC_REFERENCE = "packaged:POPC"
POPE_REFERENCE = "packaged:POPE"
POPS_REFERENCE = "packaged:POPS"
CHL_REFERENCE = "packaged:CHL"
# -------------------------------------------------------------------
# Mixed membrane construction
# -------------------------------------------------------------------
CHOLESTEROL_TRIALS_PER_SITE = 160
CHOLESTEROL_GLOBAL_ATTEMPTS = 12
CHOLESTEROL_MINIMUM_HH_DISTANCE_ANGSTROM = 1.15
MODULAR_LIPID_JUNCTION_MINIMUM_ANGSTROM = 1.25
MODULAR_LIPID_JUNCTION_MAXIMUM_ANGSTROM = 1.75
MODULAR_INTERNAL_CLEARANCE_MINIMUM_ANGSTROM = 1.75
MODULAR_INTERNAL_CLEARANCE_PREFERRED_ANGSTROM = 2.30
# -------------------------------------------------------------------
# Protein-membrane construction
# -------------------------------------------------------------------
PROTEIN_EXCLUSION_SHELL_ANGSTROM = 0.0
FOOTPRINT_Z_MARGIN_ANGSTROM = 8.0
EDGE_PADDING_ANGSTROM = 20.0
PROTEIN_MEMBRANE_INITIAL_BOX_Z_ANGSTROM = 200.0
# -------------------------------------------------------------------
# Membrane-only construction
# -------------------------------------------------------------------
MEMBRANE_ONLY_INITIAL_BOX_Z_ANGSTROM = 96.0
# -------------------------------------------------------------------
# Steric packing
# -------------------------------------------------------------------
FIXED_HEAVY_LIPID_HEAVY_ANGSTROM = 2.0
FIXED_HEAVY_LIPID_HYDROGEN_ANGSTROM = 1.2
LIPID_HYDROGEN_HYDROGEN_ANGSTROM = 1.15
# -------------------------------------------------------------------
# Ion placement
# -------------------------------------------------------------------
MINIMUM_ION_FIXED_HEAVY_DISTANCE_ANGSTROM = 4.0
MINIMUM_ION_ION_DISTANCE_ANGSTROM = 5.0
# Constants required by the active membrane builder
MODULAR_JUNCTION_MAXIMUM_ANGSTROM = 1.75
MODULAR_JUNCTION_MINIMUM_ANGSTROM = 1.25
