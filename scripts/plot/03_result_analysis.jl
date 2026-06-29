# Result analysis and plot generation script.
#
# Prerequisites: run 01_a_setup.jl → 01_b_data_split.jl → 02_location_calibration.jl
# to produce the calibrated parameter file referenced by OUT_DIR / RESULT_FN.
#
# Run in a fresh Julia session (without 01_a_setup.jl first).

using CoralBloxCalib
import CoralBloxCalib.viz
using CairoMakie: save

if @isdefined(LocationDataStore) && LocationDataStore !== CoralBloxCalib.LocationDataStore
    error(
        "LocationDataStore is already defined in this session from a different source. " *
        "Please restart Julia and run this script in a fresh session (do not run " *
        "01_a_setup.jl first)."
    )
end

_src_dir = joinpath(dirname(dirname(dirname(@__FILE__))), "src")

include(joinpath(_src_dir, "01_a_setup.jl"))
include(joinpath(_src_dir, "01_b_data_split.jl"))

# ----- Load cover, run model --------------------------------------------------

init_cover = deserialize(INIT_COVER_PATH)
construct_cover!(dom, init_cover, location_classification.consecutive_classification)

calibrated_params = deserialize(joinpath(OUT_DIR, RESULT_FN))

include(joinpath(_src_dir, "common", "params_extraction.jl"))
dom, scen = setup_run(dom, calibrated_params)
rs_raw = ADRIA.run_model(dom, scen[1, :])
mkpath(OUT_DIR)

# ----- Generate plots ---------------------------------------------------------

# Observation locations map
f_obs_loc_map = viz.plot_observation_locs(CALIBRATION_STORE, VALIDATION_STORE)
save(joinpath(OUT_DIR, "obs_loc_map.png"), f_obs_loc_map)

# Functional group cover proportions
f_taxa_cover = viz.taxa_cover_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_cov.png"), f_taxa_cover)

f_taxa_pop = viz.taxa_population_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_pop.png"), f_taxa_pop)

f_size_class = viz.temporal_size_class_proportions(rs_raw.raw; fig_size=(900, 600))
save(joinpath(OUT_DIR, "locs_size.png"), f_size_class)

# Summary stats
n_validation_locs = length(VALIDATION_STORE.ltmp_unique_ids)
validation_sccs = [
    viz.collect_error_stats(rs_raw.raw, id, dom; observations=VALIDATION_STORE).srcc
    for id in 1:n_validation_locs
]
validation_sccs_sortperm = sortperm(validation_sccs)
@info "Three highest SCC validation reefs: $(VALIDATION_STORE.ltmp_unique_ids[validation_sccs_sortperm][end-2:end])"
@info "Three lowest  SCC validation reefs: $(VALIDATION_STORE.ltmp_unique_ids[validation_sccs_sortperm][1:3])"

stats_calib = viz.collect_error_stats(rs_raw.raw, dom; observations=CALIBRATION_STORE)
stats_valid = viz.collect_error_stats(rs_raw.raw, dom; observations=VALIDATION_STORE)
n_calib = length(CALIBRATION_STORE.ltmp_unique_ids)
n_valid = length(VALIDATION_STORE.ltmp_unique_ids)
@info "Calibration locations where model outperforms benchmark: $(sum(stats_calib.rmse_model .< stats_calib.rmse_benchmark)) / $n_calib"
@info "Validation  locations where model outperforms benchmark: $(sum(stats_valid.rmse_model .< stats_valid.rmse_benchmark)) / $n_valid"
@info "Mean model calib. RMSE: $(mean(stats_calib.rmse_model))"
@info "Mean model valid.  RMSE: $(mean(stats_valid.rmse_model))"
@info "Mean model calib. SRCC: $(mean(stats_calib.srcc))"
@info "Mean model valid.  SRCC: $(mean(stats_valid.srcc))"

# Regional comparison plots (03_b)
viz.save_regional_analysis_plots(
    rs_raw.raw, dom, CALIBRATION_STORE, VALIDATION_STORE, OUT_DIR,
    ltmp_north, ltmp_central, ltmp_south,
    NORTH_MASK, CENTRAL_MASK, SOUTH_MASK
)

# Metric analysis plots (03_c)
viz.save_metric_analysis_plots(
    rs_raw.raw, dom, CALIBRATION_STORE, VALIDATION_STORE, OUT_DIR
)

# Per-location time-series plots (03_d)
cyc_scens = dom.cyclone_mortality_scens[scenarios=1]
dhw_scens  = dom.dhw_scens[scenarios=1]
disturbances = open_dataset(
    joinpath(root_path, "datasets", "ltmp_data", "disturbances.nc")
).layer

viz.save_location_timeseries_plots(
    rs_raw.raw, dom, CALIBRATION_STORE, VALIDATION_STORE, OUT_DIR,
    dhw_scens, cyc_scens, disturbances
)
