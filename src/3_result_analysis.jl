using ADRIA: RESULTS
using ProgressMeter

include("./common/common.jl")
include("./common/cover_construction.jl")
include("./plot/plot.jl")
include("1_setup.jl")
include("common/param_bounds.jl")

include("3a_analysis_setup.jl")

include("3b_regional_analysis.jl")

mkpath(OUT_DIR)

f_taxa_cover = taxa_cover_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_cov.png"), f_taxa_cover)

f_taxa_pop = taxa_population_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_pop.png"), f_taxa_pop)

f_size_class = temporal_size_class_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_size.png"), f_size_class)

mkpath("$(OUT_DIR)/loc_plots")

calibration_save_dir::String = joinpath(OUT_DIR, "calibration_locations")
validation_save_dir::String = joinpath(OUT_DIR, "validation_locations")

mkpath(calibration_save_dir)
mkpath(validation_save_dir)

# Plot Metrics
metrics_save_dir::String = joinpath(OUT_DIR, "metrics")
mkpath(metrics_save_dir)

rmse_diff_map = plot_rmse_diff_map(
    rs_raw.raw;
    observations=VALIDATION_STORE,
    fig_opts=Dict(:title => "Benchmark RMSE - Model RMSE\nValidation Locations")
)
save(joinpath(metrics_save_dir, "rmse_diff_map.png"), rmse_diff_map)

pearson_coeff_map = plot_pcc_map(
    rs_raw.raw;
    observations=VALIDATION_STORE,
    fig_opts=Dict(:title => "Pearson Correlation Coefficient\nValidation Locations")
)
save(joinpath(metrics_save_dir, "pcc_map.png"), pearson_coeff_map)

# cyc_scens = _mortality_to_cyc_category(copy(dom.cyclone_mortality_scens[scenarios=1, species=5]))
cyc_scens = dom.cyclone_mortality_scens[scenarios=1]
dhw_scens = dom.dhw_scens[scenarios=1]
disturbances_path = "C:/Users/pribeiro/AIMS/Code/ltmp_calibration/datasets/ltmp_data/disturbances.nc"
disturbances = open_dataset(disturbances_path).layer

# include("./plot/plot.jl")
f_obs_loc_map = plot_observation_locs(CALIBRATION_STORE, VALIDATION_STORE)
save(joinpath(OUT_DIR, "obs_loc_map.png"), f_obs_loc_map)

# validation_ids = ["16015100104", "23048100104", "19209100104", "14137100104"]
# loc_comparison_fig = plot_location_comparison_highlights(rs_raw.raw, validation_ids,
#     cyc_scens, dhw_scens, disturbances; observations=VALIDATION_STORE)
# save(joinpath(OUT_DIR, "loc_comparison.png"), loc_comparison_fig)

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
        observations=VALIDATION_STORE, hide_disturbances=true)
    save(joinpath(validation_save_dir, "loc_$(reef_id).png"), f)
end

rmse_diff_validation = sort(rmse_diff(rs_raw.raw, VALIDATION_STORE))
f_rmse_diff_validation = plot_rmse_scatter(rmse_diff_validation, "Validation")
save(joinpath(metrics_save_dir, "rmse_diff_validation.png"), f_rmse_diff_validation)

rmse_diff_calibration = sort(rmse_diff(rs_raw.raw, CALIBRATION_STORE))
f_rmse_diff_calibration = plot_rmse_scatter(rmse_diff_calibration, "Calibration")
save(joinpath(metrics_save_dir, "rmse_diff_calibration.png"), f_rmse_diff_calibration)

pcc_validation = pcc_locs(rs_raw.raw, VALIDATION_STORE)
f_pcc_validation = plot_pcc_scatter(pcc_validation)
save(joinpath(metrics_save_dir, "pcc_validation.png"), f_pcc_validation)

pcc_calibration = pcc_locs(rs_raw.raw, CALIBRATION_STORE)
f_pcc_calibration = plot_pcc_scatter(pcc_calibration)
save(joinpath(metrics_save_dir, "pcc_calibration.png"), f_pcc_calibration)

# ! function out

raw_data = rs_raw.raw
n_validation_obs = length(VALIDATION_STORE.ltmp_unique_ids)
ltmp_loc_indexes = collect(1:n_validation_obs)
error_stats = collect_error_stats.(Ref(raw_data), ltmp_loc_indexes; observations=VALIDATION_STORE)
rmse_, benchmark_, cc_, maee_, bias_ = eachrow(hcat(map(collect, error_stats)...))

x_labels = ["Model RMSE", "Benchmark RMSE", "PCC", "MAE", "Bias"]
x = 1:length(x_labels)
validation_ids = [findfirst(id .== VALIDATION_STORE.domain_gpkg.UNIQUE_ID) for id in VALIDATION_STORE.ltmp_unique_ids]
y_coords = VALIDATION_STORE.domain_gpkg[validation_ids, "Y_COORD"]
y_coords_sortperm = sortperm(y_coords, rev=false)
y_coords_sortperm_rev = sortperm(y_coords, rev=true)
fig = Figure()
ax = Axis(
    fig[1, 1],
    title="Calibration Metrics vs Latitude",
    xlabel="",
    ylabel="Latitude",
    xticks=(x, x_labels),
    yticks=(1:length(y_coords), string.(y_coords[y_coords_sortperm])),
    xticklabelrotation=(π / 6)
)
heat = hcat(rmse_, benchmark_, cc_, maee_, bias_)[y_coords_sortperm_rev, :]
hm = heatmap!(ax, x, 1:length(y_coords), heat')
Colorbar(fig[:, end+1], hm)
fig

# !

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

rmse_ = 0.0
benchmark_ = 0.0
cc_ = 0.0
maee_ = 0.0
bias_ = 0.0

@showprogress desc = "Calculating calibration error." for i in 1:length(CALIBRATION_STORE.ltmp_cover_to_domain)
    local tmp = collect_error_stats(rs_raw.raw, i; observations=CALIBRATION_STORE)
    global rmse_ += tmp[1] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
    global benchmark_ += tmp[2] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
    global cc_ += tmp[3] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
    global maee_ += tmp[4] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
    global bias_ += tmp[5] / length(CALIBRATION_STORE.ltmp_cover_to_domain)
end

rmse_ = 0.0
benchmark_ = 0.0
cc_ = 0.0
maee_ = 0.0
bias_ = 0.0
count_ = 0

@showprogress desc = "Calculating validation error." for i in 1:length(VALIDATION_STORE.ltmp_cover_to_domain)
    local tmp = collect_error_stats(rs_raw.raw, i; observations=VALIDATION_STORE)
    if tmp[1] < tmp[2]
        global count_ += 1
    end
    global rmse_ += tmp[1] / length(VALIDATION_STORE.ltmp_cover_to_domain)
    global benchmark_ += tmp[2] / length(VALIDATION_STORE.ltmp_cover_to_domain)
    global cc_ += tmp[3] / length(VALIDATION_STORE.ltmp_cover_to_domain)
    global maee_ += tmp[4] / length(VALIDATION_STORE.ltmp_cover_to_domain)
    global bias_ += tmp[5] / length(VALIDATION_STORE.ltmp_cover_to_domain)
end

mean(getindex.(collect_error_stats.([rs_raw.raw], collect(1:length(VALIDATION_STORE.ltmp_cover_to_domain)); observations=VALIDATION_STORE), 1))
mean(getindex.(collect_error_stats.([rs_raw.raw], collect(1:length(VALIDATION_STORE.ltmp_cover_to_domain)); observations=VALIDATION_STORE), 2))
