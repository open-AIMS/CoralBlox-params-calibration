using CairoMakie
using GeoMakie

function plot_rmse_scatter(rmse_diff, observation_type::String)
    _rmse_diff = sort(rmse_diff)

    # Calculate mean and confidence intervals
    mean_rmse_diff = mean(_rmse_diff)
    std_rmse_diff = std(_rmse_diff)
    # ci_lower = mean_rmse_diff - 1.96 * std_rmse_diff / sqrt(length(_rmse_diff))
    # ci_upper = mean_rmse_diff + 1.96 * std_rmse_diff / sqrt(length(_rmse_diff))

    n_greater_than_zero = sum(_rmse_diff .> 0)
    success_rate = round((n_greater_than_zero / length(_rmse_diff)) * 100, digits=2)

    # Create scatter plot
    fig = Figure()
    ax = Axis(fig[1, 1], title="(Benchmark RMSE - Model RMSE) for $(titlecase(observation_type)) Locations\n
    $(success_rate)% of positive diffs", xlabel="Index", ylabel="RMSE Difference")
    scatter!(ax, 1:length(_rmse_diff), _rmse_diff, color=:blue, label="RMSE Diff")

    # Add horizontal line for mean
    hlines!(ax, [mean_rmse_diff], color=:red, linestyle=:dash, label="Mean")

    # Add confidence interval bands
    #band!(ax, collect(1:length(_rmse_diff)),
    #    fill(ci_lower, length(_rmse_diff)), fill(ci_upper, length(_rmse_diff)),
    #    color=(1, 0, 0, 0.8), label="95% CI")

    axislegend(ax)
    fig
end

function plot_pcc_scatter(pcc)
    lower_bound, mean, upper_bound = average_cc(pcc)

    # Create scatter plot
    fig = Figure()

    _yticks = sort(unique(vcat((-1:0.5:1), round(mean, digits=2))))
    ax = Axis(
        fig[1, 1],
        title="Pearson Correlation Coefficient (PCC) Scatter Plot",
        xlabel="Index",
        ylabel="PCC",
        yticks=(_yticks)
    )

    scatter!(ax, 1:length(pcc), sort(pcc), color=:blue, label="PCC")

    # Add horizontal line for mean
    hlines!(ax, [mean], color=:red, linestyle=:dash, label="Mean")

    # Add confidence interval bands
    # band!(ax, collect(1:length(pcc)),
    #     fill(lower_bound, length(pcc)), fill(upper_bound, length(pcc)),
    #     color=(0, 0, 1, 0.2), label="Bounds")

    legend_els = [
        MarkerElement(color=:blue, marker=:circle),
        LineElement(color=:red, linestyle=:dash),
    ]

    Legend(
        fig[1, 2],
        legend_els,
        ["PCC", "Mean"]
    )

    fig
end

function plot_rmse_diff_map(
    raw_data::Array{Float64,3};
    fig_opts::Dict=Dict(),
    observations::LocationDataStore=COMBINED_STORE,
)::Figure
    n_validation_obs = length(observations.ltmp_unique_ids)
    ltmp_loc_indexes = collect(1:n_validation_obs)

    error_stats = collect_error_stats.([raw_data], ltmp_loc_indexes; observations=observations)
    rmse_, benchmark_, cc_, maee_, bias_ = eachrow(hcat(map(collect, error_stats)...))
    rmse_diffs = benchmark_ .- rmse_

    fig_size = get(fig_opts, :size, FIG_SIZE[:map])
    fig_title = get(fig_opts, :title, "Benchmark RMSE - Model RMSE")
    fig = Figure(size=fig_size)
    ax = GeoAxis(
        fig[1, 1],
        dest="+proj=latlong +datum=WGS84",
        title=fig_title,
        titlesize=18,
        xlabelsize=15,
        ylabelsize=15,
        xlabel="Longitude",
        ylabel="Latitude"
    )
    domain_gpkg = observations.domain_gpkg
    observation_gpkg = domain_gpkg[domain_gpkg.UNIQUE_ID.∈[observations.ltmp_unique_ids], :]

    poly!.(ax, domain_gpkg.geom, color=:gray)
    poly!.(ax, observation_gpkg.geom, color=:red)
    scatter!(ax, observation_gpkg.X_COORD, observation_gpkg.Y_COORD;
        markersize=35, color=rmse_diffs, colormap=:bam, alpha=0.8,
        strokewidth=1, strokecolor=(:gray, 0.1)
    )
    max_val, min_val = extrema(rmse_diffs)
    up_limit = maximum(abs.((max_val, min_val)))
    lower_limit = -up_limit
    Colorbar(
        fig[1, 2];
        limits=(lower_limit, up_limit), colormap=:bam, label="Model RMSE - Benchmark RMSE"
    )

    return fig
end

function plot_pcc_map(
    raw_data::Array{Float64,3};
    fig_opts::Dict=Dict(),
    observations::LocationDataStore=COMBINED_STORE,
)::Figure
    n_validation_obs = length(observations.ltmp_unique_ids)
    ltmp_loc_indexes = collect(1:n_validation_obs)

    error_stats = collect_error_stats.(Ref(raw_data), ltmp_loc_indexes; observations=observations)
    rmse_, benchmark_, cc_, maee_, bias_ = eachrow(hcat(map(collect, error_stats)...))

    fig_size = get(fig_opts, :size, FIG_SIZE[:map])
    fig_title = get(fig_opts, :title, "Pearson Correlation Coefficient")
    fig = Figure(size=fig_size)
    ax = GeoAxis(
        fig[1, 1],
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        title=fig_title,
        titlesize=18,
        # xlabelsize=15,
        # ylabelsize=15,
    )
    domain_gpkg = observations.domain_gpkg
    observation_gpkg = domain_gpkg[domain_gpkg.UNIQUE_ID.∈[observations.ltmp_unique_ids], :]
    poly!.(ax, domain_gpkg.geom, color=:gray)

    scatter!(ax, observation_gpkg.X_COORD, observation_gpkg.Y_COORD, markersize=35,
        color=cc_, colormap=:bam, alpha=0.8, strokewidth=1, strokecolor=(:gray, 0.1))
    max_val, min_val = extrema(cc_)
    up_limit = maximum(abs.((max_val, min_val)))
    lower_limit = -up_limit
    Colorbar(fig[1, 2], limits=(lower_limit, up_limit), colormap=:bam,
        label="Pearson Correlation Coefficient")

    return fig
end
