# Result analysis and plot generation script.
#
# Prerequisites: run scripts/run_loc_calib.jl first to produce the calibrated parameter
# file at config.out_dir/results.dat.

using ADRIA
using Statistics
using Serialization
using YAXArrays

using CoralBloxCalib
import CoralBloxCalib.viz
using CoralBloxCalib.common
using CoralBloxCalib.calibration
using CairoMakie: save

config = load_config()
dom = load_domain(config)
location_classification = load_location_classification(config.loc_class_path)
regional_data = viz.load_regional_analysis_data(
    dom, config.ltmp_modelled_obs_path, config.ltmp_shp_path
)

OUT_DIR = config.out_dir

# ----- Build parameter config and data stores ---------------------------------

cfg = CalibConfig(dom)

calib_data = build_calibration_data(
    dom,
    config.ltmp_reef_data_path, config.composition_path
)

# ----- Load cover, run model --------------------------------------------------

init_cover = deserialize(config.init_cover_path)
construct_cover!(dom, init_cover, location_classification.consecutive_classification)

calibrated_params = deserialize(joinpath(config.out_dir, "results.dat"))

params_dir_path = joinpath(OUT_DIR, "params")
mkpath(params_dir_path)
savedataset(
    build_params_dataset(
        calibrated_params, cfg.param_idxs, cfg.coral_param_names,
        cfg.biogroups_ordering, cfg.growth_accel_names
    );
    path=joinpath(params_dir_path, "calibrated_params.nc"),
    driver=:netcdf,
    overwrite=true
)
savedataset(
    build_init_cover_dataset(
        dom, init_cover, location_classification.consecutive_classification,
        calibrated_params, cfg.param_idxs, calib_data.combined_store, cfg.biogroups_ordering
    );
    path=joinpath(params_dir_path, "historic_init_cover.nc"),
    driver=:netcdf,
    overwrite=true
)

dom, scen = setup_run(
    dom,
    calibrated_params;
    param_names=cfg.coral_param_names,
    growth_accel_names=cfg.growth_accel_names,
    param_idxs=cfg.param_idxs,
    observations=calib_data.combined_store,
    biogroup_ord=cfg.biogroups_ordering,
)
rs_raw = ADRIA.run_model(dom, scen[1, :])
mkpath(OUT_DIR)

# ----- Generate plots ---------------------------------------------------------

# Observation locations map
f_obs_loc_map = viz.plot_observation_locs(calib_data.calibration_store, calib_data.validation_store)
save(joinpath(OUT_DIR, "obs_loc_map.png"), f_obs_loc_map)

# Functional group cover proportions
f_taxa_cover = viz.taxa_cover_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_cov.png"), f_taxa_cover)

f_taxa_pop = viz.taxa_population_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_pop.png"), f_taxa_pop)

f_size_class = viz.temporal_size_class_proportions(rs_raw.raw; fig_size=(900, 600))
save(joinpath(OUT_DIR, "locs_size.png"), f_size_class)

# Summary stats
n_validation_locs = length(calib_data.validation_store.ltmp_unique_ids)
validation_sccs = [
    viz.collect_error_stats(rs_raw.raw, id, dom; observations=calib_data.validation_store).srcc
    for id in 1:n_validation_locs
]
validation_sccs_sortperm = sortperm(validation_sccs)
@info "Three highest SCC validation reefs: $(calib_data.validation_store.ltmp_unique_ids[validation_sccs_sortperm][end-2:end])"
@info "Three lowest  SCC validation reefs: $(calib_data.validation_store.ltmp_unique_ids[validation_sccs_sortperm][1:3])"

stats_calib = viz.collect_error_stats(rs_raw.raw, dom; observations=calib_data.calibration_store)
stats_valid = viz.collect_error_stats(rs_raw.raw, dom; observations=calib_data.validation_store)
n_calib = length(calib_data.calibration_store.ltmp_unique_ids)
n_valid = length(calib_data.validation_store.ltmp_unique_ids)
@info "Calibration locations where model outperforms benchmark: $(sum(stats_calib.rmse_model .< stats_calib.rmse_benchmark)) / $n_calib"
@info "Validation  locations where model outperforms benchmark: $(sum(stats_valid.rmse_model .< stats_valid.rmse_benchmark)) / $n_valid"
@info "Mean model calib. RMSE: $(mean(stats_calib.rmse_model))"
@info "Mean model valid.  RMSE: $(mean(stats_valid.rmse_model))"
@info "Mean model calib. SRCC: $(mean(stats_calib.srcc))"
@info "Mean model valid.  SRCC: $(mean(stats_valid.srcc))"

# Regional comparison plots (03_b)
viz.save_regional_analysis_plots(
    rs_raw.raw, dom, calib_data.calibration_store, calib_data.validation_store, OUT_DIR,
    regional_data.ltmp_north, regional_data.ltmp_central, regional_data.ltmp_south,
    regional_data.north_mask, regional_data.central_mask, regional_data.south_mask
)

# Metric analysis plots (03_c)
viz.save_metric_analysis_plots(
    rs_raw.raw, dom, calib_data.calibration_store, calib_data.validation_store, OUT_DIR
)

# Per-location time-series plots (03_d)
root_path = dirname(dirname(dirname(@__FILE__)))
cyc_scens = dom.cyclone_mortality_scens[scenarios=1]
dhw_scens  = dom.dhw_scens[scenarios=1]
disturbances = open_dataset(
    joinpath(root_path, "datasets", "ltmp_data", "disturbances.nc")
).layer

viz.save_location_timeseries_plots(
    rs_raw.raw, dom, calib_data.calibration_store, calib_data.validation_store, OUT_DIR,
    dhw_scens, cyc_scens, disturbances
)
nothing
