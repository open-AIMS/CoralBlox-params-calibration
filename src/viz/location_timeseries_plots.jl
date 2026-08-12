"""
    save_location_timeseries_plots(raw_data, dom, calibration_store, test_store, out_dir,
        dhw_scens, cyc_scens, disturbances)

Generate and save per-location time-series comparison plots for all calibration and
test locations under `out_dir/calibration_locations/` and `out_dir/test_locations/`.
"""
function save_location_timeseries_plots(
    raw_data::Array{Float64,4},
    dom,
    calibration_store::LocationDataStore,
    test_store::LocationDataStore,
    out_dir::String,
    dhw_scens,
    cyc_scens,
    disturbances,
)::Nothing
    calibration_save_dir = joinpath(out_dir, "calibration_locations")
    test_save_dir = joinpath(out_dir, "test_locations")
    mkpath(calibration_save_dir)
    mkpath(test_save_dir)

    plot_type = ["calibration", "test"]
    plot_stores = [calibration_store, test_store]
    plot_save_dirs = [calibration_save_dir, test_save_dir]

    for (idx_n_locs, store) in enumerate(plot_stores)
        n_locs = length(store.ltmp_cover_to_domain)
        @showprogress desc = "Plotting $(plot_type[idx_n_locs]) locations." for i in 1:n_locs
            reef_id = store.ltmp_unique_ids[i]
            f = plot_location_comparison(
                raw_data, i, dhw_scens, cyc_scens, disturbances;
                dom=dom,
                fig_opts=Dict{Symbol,Any}(
                    :size => (600, 600),
                    :titlesize => 10pt,
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
                observations=store
            )
            save(joinpath(plot_save_dirs[idx_n_locs], "loc_$(reef_id).png"), f)
        end
    end

    @info "Location time series plots saved to: $calibration_save_dir and $test_save_dir"
    return nothing
end
