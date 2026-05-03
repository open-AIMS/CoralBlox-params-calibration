using ADRIA: RESULTS
using ProgressMeter

include("./plot/plot.jl")
include("common/param_bounds.jl")

include("03_a_analysis_setup.jl")
include("03_b_regional_analysis.jl")
include("03_c_metric_analysis.jl")
include("03_d_locs_time_series.jl")

# Map with observation and validation locations
f_obs_loc_map = plot_observation_locs(CALIBRATION_STORE, VALIDATION_STORE)
save(joinpath(OUT_DIR, "obs_loc_map.png"), f_obs_loc_map)

# Functional group cover proportions
f_taxa_cover = taxa_cover_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_cov.png"), f_taxa_cover)

# Functional group population proportions
f_taxa_pop = taxa_population_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_pop.png"), f_taxa_pop)

# Multiplot (one by functional group): population proportion per size class
f_size_class = temporal_size_class_proportions(rs_raw.raw; fig_size=(900, 600))
save(joinpath(OUT_DIR, "locs_size.png"), f_size_class)

# Identify 3 higher and 3 lower scc locations
n_validation_locs = length(VALIDATION_STORE.ltmp_unique_ids)
validation_sccs = [collect_error_stats(rs_raw.raw, id; observations=VALIDATION_STORE)[6] for id in 1:n_validation_locs]
validation_sccs_sortperm = sortperm(validation_sccs)
validation_worst_3_scc = validation_sccs[validation_sccs_sortperm][1:3]
validation_worst_3_scc_idx = VALIDATION_STORE.ltmp_unique_ids[validation_sccs_sortperm][1:3]
validation_best_3_scc = validation_sccs[validation_sccs_sortperm][end-2:end]
validation_best_3_scc_idx = VALIDATION_STORE.ltmp_unique_ids[validation_sccs_sortperm][end-2:end]

@info "Three highest SCC validation reefs are: $validation_best_3_scc_idx"
@info "Three lowest SCC validation reefs are: $validation_worst_3_scc_idx"

# Print stats
stats_calib = collect_error_stats(rs_raw.raw; observations=CALIBRATION_STORE)
n_calib_locs = length(CALIBRATION_STORE.ltmp_unique_ids)
n_model_outperformed_calib = sum(stats_calib.rmse_model .< stats_calib.rmse_benchmark)

stats_valid = collect_error_stats(rs_raw.raw; observations=VALIDATION_STORE)
n_valid_locs = length(VALIDATION_STORE.ltmp_unique_ids)
n_model_outperformed_valid = sum(stats_valid.rmse_model .< stats_valid.rmse_benchmark)

n_locs = n_valid_locs + n_calib_locs
n_model_outperformed = n_model_outperformed_valid + n_model_outperformed_calib
@info "Number of calibration locations where model outperformed benchmark: $(n_model_outperformed_calib) of $(n_calib_locs)"
@info "Number of validation locations where model outperformed benchmark: $(n_model_outperformed_valid) of $(n_valid_locs)"
@info "Total number of locations where model outperformed benchmark: $(n_model_outperformed) of $(n_locs)"

# Error stats for location with UNIQUE_ID = uid_
# uid_ = "20348100104"
# collect_error_stats(rs_raw.raw, findfirst(VALIDATION_STORE.ltmp_unique_ids .== uid_); observations=VALIDATION_STORE)

# Mean model metrics
calibration_error_stats = collect_error_stats(
    rs_raw.raw; observations=CALIBRATION_STORE
)
validation_error_stats = collect_error_stats(
    rs_raw.raw; observations=VALIDATION_STORE
)

@info "Mean model calib.RMSE: $(mean(calibration_error_stats.rmse_model))"
@info "Mean model valid.  RMSE: $(mean(validation_error_stats.rmse_model))"

@info "Mean model calib. SRCC: $(mean(calibration_error_stats.srcc))"
@info "Mean model valid. SRCC: $(mean(validation_error_stats.srcc))"
