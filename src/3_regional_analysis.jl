using ADRIA: RESULTS
using ProgressMeter

include("./common/common.jl")
include("./common/cover_construction.jl")
include("./common/plotting/plotting.jl")
include("1_setup.jl")
include("common/param_bounds.jl")

init_cover = deserialize(INIT_COVER_PATH)

construct_cover!(
    dom, init_cover, location_classification.consecutive_classification
)

# ----- LOAD CALIBRATED RESULTS -----

coral_param_fn = joinpath(OUT_DIR, RESULT_FN)
coral_params = deserialize(coral_param_fn)

# Load values into scenario dataframe

dom, scen = setup_run(
    dom,
    coral_params
)

rs_raw = ADRIA.run_model(dom, scen[1, :])

f = plot_all_regions(
    dom, rs_raw
)

save_dir = OUT_DIR

mkpath(save_dir)

save(joinpath(save_dir, "locs_reg.png"), f)

f = taxa_cover_proportions(rs_raw.raw)

save(joinpath(save_dir, "locs_taxa_cov.png"), f)

f = taxa_population_proportions(rs_raw.raw)

save(joinpath(save_dir, "locs_taxa_pop.png"), f)

f = temporal_size_class_proportions(rs_raw.raw)

save(joinpath(save_dir, "locs_size.png"), f)

mkpath("$(save_dir)/loc_plots")

calibration_save_dir::String = joinpath(OUT_DIR, "calibration_locations")
validation_save_dir::String = joinpath(OUT_DIR, "validation_locations")

mkpath(calibration_save_dir)
mkpath(validation_save_dir)

@showprogress desc="Plotting calibration locations." for i in 1:length(CALIBRATION_STORE.ltmp_cover_to_domain)
    location_comparison(rs_raw.raw, i, calibration_save_dir; observations=CALIBRATION_STORE)
end
@showprogress desc="Plotting validation locations." for i in 1:length(VALIDATION_STORE.ltmp_cover_to_domain)
    location_comparison(rs_raw.raw, i, validation_save_dir; observations=VALIDATION_STORE)
end

rmse_      = 0.0
benchmark_ = 0.0
cc_        = 0.0
maee_      = 0.0
bias_      = 0.0

@showprogress desc="Calculating calibration error." for i in 1:length(CALIBRATION_STORE.ltmp_cover_to_domain)
    local tmp = collect_error_stats(rs_raw.raw, i; observations=CALIBRATION_STORE)
    global rmse_      += tmp[1] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
    global benchmark_ += tmp[2] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
    global cc_        += tmp[3] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
    global maee_      += tmp[4] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
    global bias_      += tmp[5] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
end

rmse_      = 0.0
benchmark_ = 0.0
cc_        = 0.0
maee_      = 0.0
bias_      = 0.0
count_     = 0

@showprogress desc="Calculating validation error." for i in 1:length(VALIDATION_STORE.ltmp_cover_to_domain)
    local tmp = collect_error_stats(rs_raw.raw, i; observations=VALIDATION_STORE)
    if tmp[1] < tmp[2]
    global count_ += 1
end
    global rmse_      += tmp[1] / length(VALIDATION_STORE.ltmp_cover_to_domain)
    global benchmark_ += tmp[2] / length(VALIDATION_STORE.ltmp_cover_to_domain)
    global cc_        += tmp[3] / length(VALIDATION_STORE.ltmp_cover_to_domain)
    global maee_      += tmp[4] / length(VALIDATION_STORE.ltmp_cover_to_domain)
    global bias_      += tmp[5] / length(VALIDATION_STORE.ltmp_cover_to_domain)
end
