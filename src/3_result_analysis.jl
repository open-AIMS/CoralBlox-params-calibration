using ADRIA: RESULTS
using ProgressMeter

include("./common/common.jl")
include("./common/cover_construction.jl")
include("./plot/plot.jl")
include("1_setup.jl")
include("common/param_bounds.jl")

include("3a_analysis_setup.jl")

include("3b_regional_analysis.jl")

include("3c_metric_analysis.jl")

f_taxa_cover = taxa_cover_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_cov.png"), f_taxa_cover)

f_taxa_pop = taxa_population_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_pop.png"), f_taxa_pop)

f_size_class = temporal_size_class_proportions(rs_raw.raw; fig_size=(900, 600))
save(joinpath(OUT_DIR, "locs_size.png"), f_size_class)

calibration_save_dir::String = joinpath(OUT_DIR, "calibration_locations")
validation_save_dir::String = joinpath(OUT_DIR, "validation_locations")

mkpath(calibration_save_dir)
mkpath(validation_save_dir)


cyc_scens = dom.cyclone_mortality_scens[scenarios=1]
dhw_scens = dom.dhw_scens[scenarios=1]
disturbances_path = "C:/Users/pribeiro/AIMS/Code/ltmp_calibration/datasets/ltmp_data/disturbances.nc"
disturbances = open_dataset(disturbances_path).layer

f_obs_loc_map = plot_observation_locs(CALIBRATION_STORE, VALIDATION_STORE)
save(joinpath(OUT_DIR, "obs_loc_map.png"), f_obs_loc_map)

n_calibration_locs = length(CALIBRATION_STORE.ltmp_cover_to_domain)
@showprogress desc = "Plotting calibration locations." for i in 1:n_calibration_locs
    reef_id = CALIBRATION_STORE.ltmp_unique_ids[i]
    f = plot_location_comparison(rs_raw.raw, i, dhw_scens, cyc_scens, disturbances;
        observations=CALIBRATION_STORE)
    save(joinpath(calibration_save_dir, "loc_$(reef_id).png"), f)
end

n_validation_locs = length(VALIDATION_STORE.ltmp_cover_to_domain)
@showprogress desc = "Plotting validation locations." for i in 1:n_validation_locs
    reef_id = VALIDATION_STORE.ltmp_unique_ids[i]
    f = plot_location_comparison(rs_raw.raw, i, dhw_scens, cyc_scens, disturbances;
        observations=VALIDATION_STORE)
    save(joinpath(validation_save_dir, "loc_$(reef_id).png"), f)
end



# * 3 best and 3 worst
# ! Identify 3 best and 3 worst locations
n_validation_locs = length(VALIDATION_STORE.ltmp_unique_ids)
# validaiton_ids = [findfirst(id .== VALIDATION_STORE.domain_gpkg.UNIQUE_ID) for id in VALIDATION_STORE.ltmp_unique_ids]
validation_pccs = [collect_error_stats(rs_raw.raw, id; observations=VALIDATION_STORE)[3] for id in 1:n_validation_locs]
validation_pccs_sortperm = sortperm(validation_pccs)
validation_worst_3_pcc = validation_pccs[validation_pccs_sortperm][1:3]
validation_worst_3_pcc_idx = VALIDATION_STORE.ltmp_unique_ids[validation_pccs_sortperm][1:3]
validation_best_3_pcc = validation_pccs[validation_pccs_sortperm][end-2:end]
validation_best_3_pcc_idx = VALIDATION_STORE.ltmp_unique_ids[validation_pccs_sortperm][end-2:end]
# ! -------------------------------------

# rmse_ = 0.0
# benchmark_ = 0.0
# cc_ = 0.0
# maee_ = 0.0
# bias_ = 0.0

# @showprogress desc = "Calculating calibration error." for i in 1:length(CALIBRATION_STORE.ltmp_cover_to_domain)
#     local tmp = collect_error_stats(rs_raw.raw, i; observations=CALIBRATION_STORE)
#     global rmse_ += tmp[1] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
#     global benchmark_ += tmp[2] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
#     global cc_ += tmp[3] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
#     global maee_ += tmp[4] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
#     global bias_ += tmp[5] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
# end

# rmse_ = 0.0
# benchmark_ = 0.0
# cc_ = 0.0
# maee_ = 0.0
# bias_ = 0.0
# count_ = 0

# @showprogress desc = "Calculating validation error." for i in 1:length(VALIDATION_STORE.ltmp_cover_to_domain)
#     local tmp = collect_error_stats(rs_raw.raw, i; observations=VALIDATION_STORE)
#     if tmp[1] < tmp[2]
#         global count_ += 1
#     end
#     global rmse_ += tmp[1] / length(VALIDATION_STORE.ltmp_cover_to_domain)
#     global benchmark_ += tmp[2] / length(VALIDATION_STORE.ltmp_cover_to_domain)
#     global cc_ += tmp[3] / length(VALIDATION_STORE.ltmp_cover_to_domain)
#     global maee_ += tmp[4] / length(VALIDATION_STORE.ltmp_cover_to_domain)
#     global bias_ += tmp[5] / length(VALIDATION_STORE.ltmp_cover_to_domain)
# end
