calibration_save_dir::String = joinpath(OUT_DIR, "calibration_locations")
validation_save_dir::String = joinpath(OUT_DIR, "validation_locations")

mkpath(calibration_save_dir)
mkpath(validation_save_dir)

cyc_scens = dom.cyclone_mortality_scens[scenarios=1]
dhw_scens = dom.dhw_scens[scenarios=1]
disturbances_path = joinpath(root_path, "datasets/ltmp_data/disturbances.nc")
disturbances = open_dataset(disturbances_path).layer

n_calibration_locs = length(CALIBRATION_STORE.ltmp_cover_to_domain)
n_validation_locs = length(VALIDATION_STORE.ltmp_cover_to_domain)

plot_type = ["calibration", "validation"]
plot_stores = [CALIBRATION_STORE, VALIDATION_STORE]
plot_save_dirs = [calibration_save_dir, validation_save_dir]
for (idx_n_locs, n_locs) in enumerate([n_calibration_locs, n_validation_locs])
    @showprogress desc = "Plotting $(plot_type[idx_n_locs]) locations." for i in 1:n_locs
        reef_id = plot_stores[idx_n_locs].ltmp_unique_ids[i]
        f = plot_location_comparison(
            rs_raw.raw, i, dhw_scens, cyc_scens, disturbances;
            fig_opts=Dict{Symbol,Any}(
                :size => (600, 600),
                :titlesize => 10pt,
                # :titlehalign => :left,
                :titlevalign => :bottom,
                :titlefont => :bold
            ),
            axis_opts=Dict{Symbol,Any}(
                :ylabelsize => 9pt,
                :xlabelsize => 9pt,
                :xticklabelsize => 9pt,
                :yticklabelsize => 9pt,
                :model_dist_yticks => (0:2:10),
                :model_dist_limits => (nothing, (0.0, 10)),
                :benthic_yticks => (0:0.2:1, string.(0:20:100) .* "%"),
                :benthic_limits => (nothing, (0.0, 0.6)),
                :benthic_ylabel => "Coral composition",
                :xticks => (2008:2022, string.(2008:2022)),
                :xticklabelrotation => π / 3,
                :model_vs_obs_yticks => (0:0.2:1, string.(0:20:100) .* "%"),
                :model_vs_obs_limits => (nothing, (0.0, 0.5)),
                :model_vs_obs_ylabel => "Relative coral cover",
            ),
            opts=Dict{Symbol,Any}(
                :model_vs_obs_markersize => 10,
                :model_vs_obs_linewidth => 2.5,
                :benthic_linewidth => 2.5,
                :model_dist_linewidth => 2.5,
                :show_ltmp_dist => false,
            ),
            observations=plot_stores[idx_n_locs]
        )
        save(joinpath(plot_save_dirs[idx_n_locs], "loc_$(reef_id).png"), f)
    end
end

@info "Location time series plots saved to: $calibration_save_dir and $validation_save_dir"
