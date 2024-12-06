using CairoMakie

function diff_score_barplot(
    score_a::Vector{Float64},
    score_b::Vector{Float64},
    n_observations::Vector{Int64};
    fig_size::NTuple{2,Int64}=(1000, 800),
    title::String="",
    xlabel::String="",
    colors::Vector{Symbol}=[:purple, :orange],
    legend_labels::Vector{String}=["", ""],
)
    score_diffs = score_a .- score_b
    n_locations::Int64 = length(n_observations)
    ylim_lower = floor(minimum([score_diffs..., n_observations...])) - 1
    ylim_upper = ceil(maximum([score_diffs..., n_observations...])) + 1
    limits = (nothing, (ylim_lower, ylim_upper))

    fig = Figure(size=fig_size)


    n_locs_1::Int64 = Int(floor(n_locations / 2))
    n_locs_2::Int64 = n_locations - n_locs_1

    n_obs_sortperm = sortperm(n_observations)
    rme_scores_sorted = score_a[n_obs_sortperm]
    benchmark_scores_sorted = score_b[n_obs_sortperm]
    n_observations_sorted = n_observations[n_obs_sortperm]

    if n_locations > 100
        # Split data into into two horizontal plots
        ax_1 = Axis(fig[1, 1], title=title, titlesize=20, limits=limits)
        ax_2 = Axis(fig[2, 1], xlabel=xlabel, xlabelsize=18, limits=limits)
        dodge_diff_plot!(
            ax_1, n_locs_1, rme_scores_sorted[1:n_locs_1], benchmark_scores_sorted[1:n_locs_1],
            n_observations_sorted[1:n_locs_1], colors; dodge_gap=0.05, gap=0.3, direction=:y
        )
        dodge_diff_plot!(
            ax_2, n_locs_2, rme_scores_sorted[n_locs_1+1:end],
            benchmark_scores_sorted[n_locs_1+1:end], n_observations_sorted[n_locs_1+1:end],
            colors; dodge_gap=0.05, gap=0.3, direction=:y
        )
        elements = [PolyElement(polycolor=colors[i]) for i in 1:length(legend_labels)]
        Legend(fig[3, 1], elements, legend_labels, orientation=:horizontal)
    else
        ax = Axis(fig[1, 1], xlabel=xlabel, xlabelsize=18, limits=limits)
        dodge_diff_plot!(
            ax, n_locations, rme_scores_sorted, benchmark_scores_sorted,
            n_observations_sorted, colors; dodge_gap=0.05, gap=0.3, direction=:y
        )
        elements = [PolyElement(polycolor=colors[i]) for i in 1:length(legend_labels)]
        Legend(fig[2, 1], elements, legend_labels, orientation=:horizontal)
    end


    return fig
end


function dodge_diff_plot!(
    ax::Axis,
    n_locations::Int64,
    rme_scores::Vector{Float64},
    benchmark_scores::Vector{Float64},
    n_observations::Vector{Int64},
    colors::Vector{Symbol};
    dodge_gap::Float64=0.03,
    gap::Float64=0.2,
    direction::Symbol=:x
)
    n_datasets = 2
    data_group_loc = reshape(repeat(1:n_locations, 1, n_datasets)', n_datasets * n_locations)
    data = reshape(hcat(rme_scores .- benchmark_scores, n_observations)', n_datasets * n_locations)
    data_group_dodge = reshape(repeat(1:n_datasets, 1, n_locations), n_datasets * n_locations)
    colors_dodge = reshape(repeat(colors, 1, n_locations), n_datasets * n_locations)

    return barplot!(
        ax,
        data_group_loc,
        data,
        dodge=data_group_dodge,
        dodge_gap=dodge_gap,
        gap=gap,
        color=colors_dodge,
        direction=direction
    )
end
