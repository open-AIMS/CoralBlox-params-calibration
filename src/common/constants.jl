"""
Central named constants for model dimensions and parameter bounds.

This file must be included before any file that references these constants
(see `common.jl` / `plot.jl` include order). Values mirror the literals
previously hardcoded across `cover_construction.jl`, `param_bounds.jl`,
`plotting.jl`, and the numbered pipeline scripts — no numeric value should be
changed here without re-deriving every dependent reshape/index.
"""

# --- Core domain dimensions -------------------------------------------------

# Number of coral functional groups / taxa (ADRIA.functional_group_names())
const N_TAXA = 5

# Number of coral size classes (ADRIA.bin_edges() / bin_widths())
const N_SIZE_CLASSES = 7

# Number of reef locations in the RME domain
const N_LOCATIONS = 3806

# Number of timesteps in a single model run output (raw_data first dimension)
const N_TIMESTEPS = 15

# --- Derived / parameter-vector dimensions ----------------------------------

# Number of flattened coral parameters per location sample (N_TAXA * N_SIZE_CLASSES)
const N_PARAMS = N_TAXA * N_SIZE_CLASSES

# Stride of a single location's sample vector in `vec_sample`:
# 1 relative cover + N_TAXA taxonomy weights + N_TAXA size-distribution lambdas
const STRIDE = 1 + 2 * N_TAXA

# Number of growth-acceleration parameters per biogroup (steepness, height, midpoint)
const N_GROWTH_ACCEL_PARAMS = 3

# --- Parameter bound constants (param_bounds.jl) ----------------------------

# Linear extension / background mortality bounds: +/-10% of the EcoRRAP mean
const LIN_EXT_LB_FACTOR = 0.9
const LIN_EXT_UB_FACTOR = 1.1
const MB_RATE_LB_FACTOR = 0.9
const MB_RATE_UB_FACTOR = 1.1

# Linear extension scale factor bounds
const LIN_EXT_SCALE_LB = 0.7
const LIN_EXT_SCALE_UB = 1.5

# Background mortality scale factor bounds
const MB_RATE_SCALE_LB = -1.0
const MB_RATE_SCALE_UB = 1.0

# Growth acceleration parameter bounds (steepness, height, midpoint)
const GROWTH_ACCEL_STEEPNESS_BOUNDS = (-20.0, -15.0)
const GROWTH_ACCEL_HEIGHT_BOUNDS = (0.0, 2.0)
const GROWTH_ACCEL_MIDPOINT_BOUNDS = (0.0, 0.3)

# Size-class distribution (lambda) bounds
const SC_DIST_LB = 0.25
const SC_DIST_UB = 2.0
