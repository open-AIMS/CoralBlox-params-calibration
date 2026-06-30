module common

using ADRIA
using DataFrames
using StatsBase
using Statistics
using DataStructures
using CSV

import ..CoralBloxCalib:
    LocationDataStore,
    ltmp_cover_idx_to_domain,
    composition_idx_to_domain,
    get_ltmp_loc_unique_id,
    extract_param_group_idx,
    extract_param_group

include("constants.jl")
include("cover_construction.jl")
include("perf_metrics.jl")

export LocationDataStore,
       ltmp_cover_idx_to_domain,
       composition_idx_to_domain,
       get_ltmp_loc_unique_id,
       extract_param_group_idx,
       extract_param_group,
       squared_expo_cdf,
       size_class_proportion,
       size_class_distribution,
       construct_location_cover!,
       compute_cover,
       construct_cover!,
       temporal_correlation,
       temporal_correlation_penalty,
       rmse,
       bias,
       MAE,
       MAEE,
       MAEE_series,
       temporal_variability,
       calib_func,
       collect_error_stats,
       rmse_diff,
       location_correlation_coefficients,
       average_cc,
       reef_taxa_error,
       reef_error,
       average_class_cover,
       START_YEAR,
       END_YEAR,
       N_TAXA,
       N_SIZE_CLASSES,
       N_LOCATIONS,
       N_TIMESTEPS,
       N_PARAMS,
       STRIDE,
       N_GROWTH_ACCEL_PARAMS,
       LIN_EXT_LB_FACTOR,
       LIN_EXT_UB_FACTOR,
       MB_RATE_LB_FACTOR,
       MB_RATE_UB_FACTOR,
       LIN_EXT_SCALE_LB,
       LIN_EXT_SCALE_UB,
       MB_RATE_SCALE_LB,
       MB_RATE_SCALE_UB,
       GROWTH_ACCEL_STEEPNESS_BOUNDS,
       GROWTH_ACCEL_HEIGHT_BOUNDS,
       GROWTH_ACCEL_MIDPOINT_BOUNDS,
       SC_DIST_LB,
       SC_DIST_UB

end
