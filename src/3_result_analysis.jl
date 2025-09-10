using ADRIA: RESULTS
using ProgressMeter

include("./plot/plot.jl")
include("common/param_bounds.jl")

include("3a_analysis_setup.jl")

include("3b_regional_analysis.jl")

include("3c_metric_analysis.jl")

include("3d_locs_time_series.jl")

f_taxa_cover = taxa_cover_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_cov.png"), f_taxa_cover)

f_taxa_pop = taxa_population_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_pop.png"), f_taxa_pop)

f_size_class = temporal_size_class_proportions(rs_raw.raw; fig_size=(900, 600))
save(joinpath(OUT_DIR, "locs_size.png"), f_size_class)

f_obs_loc_map = plot_observation_locs(CALIBRATION_STORE, VALIDATION_STORE)
save(joinpath(OUT_DIR, "obs_loc_map.png"), f_obs_loc_map)

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
