using CairoMakie

function diff_score_barplot(
    score_a::AbstractVector{Float64},
    score_b::AbstractVector{Float64},
    n_observations::AbstractVector{Int64};
    fig_size::NTuple{2,Int64}=(1000, 800),
    title::String="",
    xlabel::String="",
    colors::Vector{Symbol}=[:purple, :orange],
    legend_labels::Vector{String}=["", ""],
)
    fig = Figure(size=fig_size)

    n_obs_sortperm = sortperm(n_observations)
    sorted_score_a = score_a[n_obs_sortperm]
    sorted_score_b = score_b[n_obs_sortperm]
    n_observations_sorted = n_observations[n_obs_sortperm]

    score_diffs = sorted_score_a .- sorted_score_b
    n_locations::Int64 = length(n_observations_sorted)
    ylim_lower = floor(minimum([score_diffs..., n_observations_sorted...])) - 1
    ylim_upper = ceil(maximum([score_diffs..., n_observations_sorted...])) + 1
    limits = (nothing, (ylim_lower, ylim_upper))

    data = vcat(score_diffs', n_observations_sorted')
    if n_locations > 100
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
        elements = [PolyElement(polycolor=colors[i]) for i in 1:length(legend_labels)]
        Legend(fig[3, 1], elements, legend_labels, orientation=:horizontal)
    else
        ax = Axis(fig[1, 1], title=title, titlesize=20, xlabel=xlabel, xlabelsize=18,
            limits=limits)
        dodge_diff_plot!(
            ax, n_locations, data, colors; dodge_gap=0.05, gap=0.3, direction=:y
        )
        elements = [PolyElement(polycolor=colors[i]) for i in 1:length(legend_labels)]
        Legend(fig[2, 1], elements, legend_labels, orientation=:horizontal)
    end

    return fig
end

function _dodge_split_xticks(
    n_ticks_1::Int64, n_ticks_2::Int64; Δxtick::Int64=10
)::Tuple
    xticks_2_labels = string.((n_ticks_1+1):Δxtick:(n_ticks_1+n_ticks_2))
    return (1:Δxtick:n_ticks_1, (1:Δxtick:n_ticks_2, xticks_2_labels))
end

function dodge_diff_plot!(
    ax::Axis,
    n_locations::Int64,
    data::Matrix{Float64},
    colors::Vector{Symbol};
    dodge_gap::Float64=0.03,
    gap::Float64=0.2,
    direction::Symbol=:x
)
    n_datasets = 2
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
