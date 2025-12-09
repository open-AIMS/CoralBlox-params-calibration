single_reef_fig = begin
    plot_locs = vcat(reefs_sorted_by_srcc[end-1:end], reefs_sorted_by_srcc[1:2])
    plot_locs_idx = [findfirst(VALIDATION_STORE.ltmp_unique_ids .== loc) for loc in plot_locs]
    plot_labels = ["A", "B", "C", "D"]
    plot_positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

    # ** Stats for
    @info trunc.(validation_error_stats.rmse_model[plot_locs_idx], digits=3)
    @info trunc.(validation_error_stats.rmse_benchmark[plot_locs_idx], digits=3)
    @info trunc.(validation_error_stats.srcc[plot_locs_idx], digits=3)
    [dom.loc_data[dom.loc_data.RME_UNIQUE_ID.==plotloc, :] for plotloc in plot_locs][4]

    # ** Plot comparison
    fig = Figure(; size=(750, 750), fontsize=9pt)
    for plot_idx in 1:length(plot_locs_idx)
        # reef_idx = 1
        target_loc_data = dom.loc_data[dom.loc_data.RME_UNIQUE_ID.==plot_locs[plot_idx], :]
        title = target_loc_data.cluster_id[1]
        plot_location_comparison!(
            fig[plot_positions[plot_idx]...],
            rs_raw.raw, plot_locs[plot_idx], dhw_scens, cyc_scens, ltmp_disturbances;
            opts=Dict{Symbol,Any}(
                :show_ltmp_dist => false,
                :show_model_dist => true,
                :show_benthic_composition => true,
                :show_legends => false,
                :model_vs_obs_markersize => 10,
                :model_vs_obs_linewidth => 2.5,
                :benthic_linewidth => 2.5,
                :model_dist_linewidth => 2.5,
            ),
            axis_opts=Dict{Symbol,Any}(
                :xticklabelsvisible => plot_positions[plot_idx][1] == 2,
                :yticklabelsvisible => plot_positions[plot_idx][2] == 1,
                :xlabelvisible => plot_positions[plot_idx][1] == 2,
                :ylabelvisible => plot_positions[plot_idx][2] == 1,
                :model_vs_obs_yticks => (0:0.2:1, string.(0:20:100) .* "%"),
                :model_vs_obs_limits => (nothing, (0.0, 0.5)),
                :model_vs_obs_ylabel => "Relative coral cover",
                :model_vs_obs_xticklabelsvisible => false,
                :model_dist_xticklabelsvisible => false,
                :model_dist_yticks => (0:2:10),
                :model_dist_limits => (nothing, (0.0, 10)),
                :benthic_yticks => (0:0.2:1, string.(0:20:100) .* "%"),
                :benthic_limits => (nothing, (0.0, 0.6)),
                :benthic_ylabel => "Coral composition",
                :xticks => (2008:2022, string.(2008:2022)),
                :xticklabelrotation => π / 3,
                :ylabelsize => 9pt,
                :xlabelsize => 9pt,
                :xticklabelsize => 9pt,
                :yticklabelsize => 9pt,
                :halign => :left
            ),
            fig_opts=Dict{Symbol,Any}(
                :titlesize => 9pt,
                :title => plot_labels[plot_idx],
                :titlehalign => :left,
                :titlevalign => :bottom,
                :titlefont => :bold
            ),
            observations=VALIDATION_STORE
        )
    end
    rowsize!.(contents(fig.layout), Ref(2), Ref(Relative(2 / 10)))
    rowsize!.(contents(fig.layout), Ref(3), Ref(Relative(3 / 10)))
    fig

    model_vs_obs_legend_els = [
        MarkerElement(color=COLORS[:model_vs_obs_color_obs], marker=:circle, markersize=10),
        LineElement(color=COLORS[:model_vs_obs_color_model], linewidth=2)
    ]

    model_dist_legend_els = [
        LineElement(linewidth=3, color=COLORS[:model_dist_color_dhw], linestyle=:dot),
        LineElement(linewidth=3, color=COLORS[:model_dist_color_dist], linestyle=:dash),
    ]

    benthic_legend_els = taxa_props_legend_elements()
    taxa = String.(ADRIA.functional_group_names())
    taxa_labels = titlecase.(join.(split.(taxa, "_"), " "))

    Legend(
        fig[:, 3],
        [model_vs_obs_legend_els, model_dist_legend_els, benthic_legend_els],
        [
            ["Observed", "Model"],
            ["DHW", "Cyclone/COTS"],
            taxa_labels
        ],
        ["Coral cover", "Disturbances", "Functional groups"];
        halign=:left, valign=:center, labelsize=9pt, titlesize=9pt, labelhalign=:left,
        gridshalign=:left, titlehalign=:left
    )
    fig
end
save(fig_path * "/single_reef.png", single_reef_fig; px_per_unit=(300 / inch))