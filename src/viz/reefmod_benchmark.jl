using CairoMakie

function diff_score_barplot(
    score_model::AbstractVector{Float64},
    score_benchmark::AbstractVector{Float64},
    n_observations::AbstractVector{Int64};
    fig_size::NTuple{2,Int64}=(1000, 800),
    title::String="",
    xlabel::String="Locations",
    colors::Vector=["#D81B60", "#1E88E5", "#FFC107"],
)
    fig = Figure(size=fig_size)

    n_obs_sortperm = sortperm(n_observations)
    sorted_score_model = score_model[n_obs_sortperm]
    sorted_score_benchmark = score_benchmark[n_obs_sortperm]
    n_observations_sorted = n_observations[n_obs_sortperm]
    score_diffs = sorted_score_model .- sorted_score_benchmark
    data = vcat(score_diffs', sorted_score_model', n_observations_sorted')

    ylim_lower = floor(minimum(data)) - 1
    ylim_upper = ceil(maximum(data)) + 1
    limits = (nothing, (ylim_lower, ylim_upper))
    legend_position = [2, 1]

    # TODO: Not urgent but we could programmatically split in more plots the more data we have
    n_locations::Int64 = length(n_observations_sorted)
    n_bars = reduce(*, size(data))
    if n_bars > 200
        # Split data into into two horizontal plots
        n_locs_1::Int64 = Int(floor(n_locations / 2))
        n_locs_2::Int64 = n_locations - n_locs_1

        xticks_1, xticks_2 = _dodge_split_xticks(n_locs_1, n_locs_2)
        ax_1 = Axis(fig[1, 1], title=title, titlesize=20, limits=limits, xticks=xticks_1)
        ax_2 = Axis(fig[2, 1], xlabel=xlabel, xlabelsize=18, limits=limits, xticks=xticks_2)
        dodge_diff_plot!(
            ax_1, n_locs_1, data[:, 1:n_locs_1], colors;
            dodge_gap=0.05, gap=0.3, direction=:y
        )
        dodge_diff_plot!(
            ax_2, n_locs_2, data[:, n_locs_1+1:end], colors;
            dodge_gap=0.05, gap=0.3, direction=:y
        )
        legend_position = [3, 1]
    else
        ax = Axis(fig[1, 1], title=title, titlesize=20, xlabel=xlabel, xlabelsize=18,
            limits=limits)
        dodge_diff_plot!(
            ax, n_locations, data, colors; dodge_gap=0.05, gap=0.3, direction=:y
        )
    end

    legend_labels = ["Model - Benchmark Scores", "Model Scores", "# of observations"]
    elements = [PolyElement(polycolor=colors[i]) for i in 1:length(legend_labels)]
    Legend(fig[legend_position...], elements, legend_labels, orientation=:horizontal)

    return fig
end

function dodge_diff_plot!(
    ax::Axis,
    n_locations::Int64,
    data::Matrix{Float64},
    colors::Vector;
    dodge_gap::Float64=0.03,
    gap::Float64=0.2,
    direction::Symbol=:x
)
    n_datasets = 3
    data_group_loc = reshape(repeat(1:n_locations, 1, n_datasets)', n_datasets * n_locations)
    reshaped_data = reshape(data, reduce(*, size(data)))
    data_group_dodge = reshape(repeat(1:n_datasets, 1, n_locations), n_datasets * n_locations)
    colors_dodge = reshape(repeat(colors, 1, n_locations), n_datasets * n_locations)

    return barplot!(
        ax,
        data_group_loc,
        reshaped_data,
        dodge=data_group_dodge,
        dodge_gap=dodge_gap,
        gap=gap,
        color=colors_dodge,
        direction=direction
    )
end

function _dodge_split_xticks(
    n_ticks_1::Int64, n_ticks_2::Int64; Δxtick::Int64=10
)::Tuple
    xticks_2_labels = string.((n_ticks_1+1):Δxtick:(n_ticks_1+n_ticks_2))
    return (1:Δxtick:n_ticks_1, (1:Δxtick:n_ticks_2, xticks_2_labels))
end

function diff_map(
    rme_scores::YAXArray, benchmark_scores::YAXArray, canonical_path::String;
    fig_size=(600, 800)
)::Figure
    fig = Figure(fig_size)
    ax = Axis(fig[1, 1])

    canonical_gpkg = GDF.read(canonical_path)

    canonical_locs = canonical_gpkg.RME_UNIQUE_ID
    data_loc_idx = [findfirst(loc .== canonical_locs) for loc in locations]
    data_gpkg = canonical_gpkg[data_loc_idx, :]

    # Plot canonical reefs
    poly!(ax, canonical_gpkg.geometry)

    xy_data = Vector.(eachrow(hcat(data_gpkg.LON, data_gpkg.LAT)))
    score_diffs = (rme_scores .- benchmark_scores).data
    positive_score_diffs = score_diffs[score_diffs.>=0]
    negative_score_diffs = score_diffs[score_diffs.<0]

    map_scatter!(ax, positive_score_diffs, xy_data[score_diffs.>=0], (:orange, 0.4))
    map_scatter!(ax, negative_score_diffs, xy_data[score_diffs.<0], (:purple, 0.4))

    legend_colors = [:orange, :purple]
    legend_labels = ["Model Score >= Benchmark Score", "Benchmark Score > Model Score"]
    elements = [PolyElement(polycolor=legend_colors[i]) for i in 1:length(legend_labels)]
    Legend(fig[2, 1], elements, legend_labels, orientation=:horizontal)

    fig
end

function map_scatter!(ax, data, coord_data, color)
    min_data, max_data = extrema(data)
    Δdata = max_data - min_data
    normalized_data_diffs = data ./ Δdata
    markersizes = bin_rescale(normalized_data_diffs, collect(range(20, 70, 10)))
    for (role_id, role) in enumerate(coord_data)
        scatter!(ax, role..., markersize=markersizes[role_id], color=color)
    end
end
