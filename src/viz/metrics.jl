using CairoMakie
using GeoMakie

function plot_rmse_scatter(
    rmse_diff;
    observation_type::String="",
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    _rmse_diff = sort(rmse_diff)

    # Calculate mean and confidence intervals
    mean_rmse_diff = trunc(mean(_rmse_diff); digits=3)
    std_rmse_diff = trunc(std(_rmse_diff); digits=3)

    n_greater_than_zero = sum(_rmse_diff .> 0)
    success_rate = round((n_greater_than_zero / length(_rmse_diff)) * 100, digits=2)

    obs_title = isempty(observation_type) ? "" : "\n$(titlecase(observation_type)) data"
    axis_opts::Dict{Symbol,Any} = Dict(
        :title => "Benchmark - Model (RMSE)" * obs_title * "\n# > 0: $n_greater_than_zero | Avg: $mean_rmse_diff | Std: $std_rmse_diff",
        :ylabel => "RMSE Difference"
    )

    opts::Dict{Symbol,Any} = Dict(:metric_label => "RMSE Diff",)

    return plot_metric_scatter(
        _rmse_diff, mean_rmse_diff;
        fig_opts=fig_opts, axis_opts=axis_opts, opts=opts
    )
end

function plot_metric_scatter(
    outcomes::Vector{Float64};
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    mean_outcome = mean(outcomes)
    _yticks = sort(unique(vcat((-1:0.5:1), round(mean_outcome, digits=2))))
    axis_opts = merge(Dict(:ylabel => "Outcome", :yticks => _yticks), axis_opts)

    return plot_metric_scatter(
        sort(outcomes), mean_outcome;
        fig_opts=fig_opts, axis_opts=axis_opts, opts=opts
    )
end

function plot_metric_scatter(
    metric_data::Vector{Float64},
    agg_data::Float64;
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)
    size = get(fig_opts, :size, (600, 400))

    title = get(axis_opts, :title, "")
    xlabel = get(axis_opts, :xlabel, "Index")
    ylabel = get(axis_opts, :ylabel, "Metric value")
    yticks = get(axis_opts, :yticks, automatic)

    fig = Figure(; size=size)
    ax = Axis(
        fig[1, 1],
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        yticks=yticks
    )

    metric_color = get(opts, :metric_color, :blue)
    agg_color = get(opts, :agg_color, :red)
    metric_label = get(opts, :metric_label, "Metric")
    agg_label = get(opts, :agg_label, "Mean")

    scatter!(ax, 1:length(metric_data), metric_data, color=metric_color)
    hlines!(ax, [agg_data], color=agg_color, linestyle=:dash)

    legend_els = [
        MarkerElement(color=metric_color, marker=:circle),
        LineElement(color=agg_color, linestyle=:dash),
    ]

    Legend(fig[1, 2], legend_els, [metric_label, agg_label])

    return fig
end

function plot_rmse_diff_map(
    raw_data::Array{Float64,4},
    dom;
    observations::LocationDataStore,
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    fig_size = get!(fig_opts, :size, FIG_SIZE[:map])
    fig = Figure(; fig_opts...)
    plot_rmse_diff_map!(
        fig[1, 1], raw_data, dom;
        observations=observations, axis_opts=axis_opts, opts=opts
    )
    return fig
end
function plot_rmse_diff_map!(
    g::Union{GridPosition,GridLayout},
    raw_data::Array{Float64,4},
    dom;
    observations::LocationDataStore,
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Nothing
    n_obs = length(observations.ltmp_unique_ids)
    ltmp_loc_indexes = collect(1:n_obs)

    error_stats = collect_error_stats(raw_data, dom; observations=observations)
    rmse_diffs = error_stats.rmse_benchmark .- error_stats.rmse_model

    domain_gpkg = observations.domain_gpkg
    geometries = domain_gpkg.geometry
    lon, lat = eachcol(obs_lon_lat(observations, domain_gpkg))

    plot_metric_map!(
        g, rmse_diffs, geometries, lon, lat;
        axis_opts=axis_opts, opts=opts
    )

    return nothing
end

function obs_lon_lat(observations, domain_gpkg)::DataFrame
    ltmp_uids_obs = observations.ltmp_unique_ids
    ltmp_idx_obs = [findfirst(domain_gpkg.UNIQUE_ID .== uid) for uid in ltmp_uids_obs]
    return domain_gpkg[ltmp_idx_obs, [:LON, :LAT]]
end

function plot_correlation_map(
    raw_data::Array{Float64,4},
    dom;
    metric_type::Symbol=:scc,
    observations::LocationDataStore,
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    fig_size = pop!(fig_opts, :size, FIG_SIZE[:map])
    fig = Figure(; size=fig_size, fig_opts)
    plot_correlation_map!(
        fig[1, 1], raw_data, dom;
        metric_type=metric_type, observations=observations, axis_opts=axis_opts, opts=opts,
    )
    return fig
end
function plot_correlation_map!(
    g::Union{GridPosition,GridLayout},
    raw_data::Array{Float64,4},
    dom;
    metric_type::Symbol=:srcc,
    observations::LocationDataStore,
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Nothing
    n_validation_obs = length(observations.ltmp_unique_ids)
    ltmp_loc_indexes = collect(1:n_validation_obs)

    error_stats = collect_error_stats(raw_data, dom; observations=observations)
    cc_ = error_stats.pcc
    srcc_ = error_stats.srcc

    _metric::Vector{Float64} = if metric_type == :pcc
        cc_
    elseif metric_type == :srcc
        srcc_
    else
        error("Invalid `metric_type`. ")
    end

    domain_gpkg = observations.domain_gpkg
    geometries = domain_gpkg.geometry
    lon, lat = eachcol(obs_lon_lat(observations, domain_gpkg))

    plot_metric_map!(
        g, _metric, geometries, lon, lat;
        axis_opts=axis_opts, opts=opts
    )

    return nothing
end

function plot_metric_map!(
    g::Union{GridPosition,GridLayout},
    metric::Vector{Float64},
    geometries::Vector,
    lon_valid::Vector{Float64},
    lat_valid::Vector{Float64};
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Nothing
    axis_opts = merge(
        Dict(
            :xlabel => get(axis_opts, :xlabel, "Longitude"),
            :ylabel => get(axis_opts, :ylabel, "Latitude"),
            :title => get(axis_opts, :title, "Benchmark RMSE - Model RMSE"),
            :dest => "+proj=latlong +datum=WGS84",
        ),
        axis_opts
    )
    ax = GeoAxis(g[1, 1], ; axis_opts...)

    poly!(ax, geometries, color=:gray)


    poly!(ax, GBRMPA_MAINLAND_POLYS, color="#121212")


    max_val, min_val = extrema(metric)
    up_limit = maximum(abs.((max_val, min_val)))
    lower_limit = -up_limit

    colormap = get(opts, :colormap, :bam)
    colorrange = get(opts, :colorrange, (lower_limit, up_limit))
    alpha = get(opts, :alpha, 0.8)
    strokewidth = get(opts, :strokewidth, 1)
    strokecolor = get(opts, :strokecolor, (:gray, 0.1))
    markersize = get(opts, :markersize, 35)

    high_mask = metric .> 0
    low_mask = metric .<= 0

    scatter!(ax, lon_valid[low_mask], lat_valid[low_mask];
        colormap=colormap,
        colorrange=colorrange,
        color=metric[low_mask],
        alpha=alpha,
        marker=:circle,
        strokewidth=strokewidth,
        strokecolor=strokecolor,
        markersize=markersize,
    )

    scatter!(ax, lon_valid[high_mask], lat_valid[high_mask];
        colormap=colormap,
        colorrange=colorrange,
        color=metric[high_mask],
        alpha=alpha,
        marker=:star5,
        strokewidth=strokewidth,
        strokecolor=strokecolor,
        markersize=markersize,
    )

    if get(opts, :colorbar_visible, true)
        colorbar_label = get(opts, :colorbar_label, "Benchmark RMSE - Model RMSE")
        colorbar_ticklabelsize = get(opts, :colorbar_ticklabelsize, 16)
        colorbar_ticks = get(opts, :colorbar_ticks, automatic)
        colorbar_vertical = get(opts, :colorbar_vertical, true)
        position = colorbar_vertical ? (1, 2) : (2, 1)
        Colorbar(
            g[position...];
            colorrange=colorrange,
            colormap=colormap,
            ticks=colorbar_ticks,
            label=colorbar_label,
            ticklabelsize=colorbar_ticklabelsize,
            vertical=colorbar_vertical
        )
        # ? Can we infer the aspect ration from the LON/LAT and update the Aspect accordingly?
        # colsize!(g.layout, 1, Aspect(1, 0.15))
    end

    return nothing
end

function plot_metrics_heatmap(rs_raw, dom; fig_size=(500, 700), observations::LocationDataStore)::Figure
    error_stats = collect_error_stats(rs_raw, dom; observations=observations)
    Δrmse = error_stats.rmse_benchmark .- error_stats.rmse_model

    validation_ids = [
        findfirst(id .== observations.domain_gpkg.UNIQUE_ID)
        for id in observations.ltmp_unique_ids
    ]
    y_coords = try
        observations.domain_gpkg[validation_ids, "Y_COORD"]
    catch
        observations.domain_gpkg[validation_ids, "LAT"]
    end
    y_coords_sortperm = sortperm(y_coords, rev=false)

    heat_data = hcat(Δrmse, error_stats.pcc, error_stats.srcc, error_stats.bias)[y_coords_sortperm, :]
    x_labels = ["Benchmark - Model (RMSE)", "PCC", "SRCC", "Bias"]
    colormaps = [:Spectral, :Spectral, :Spectral, :Spectral]

    bias_up_limit = ceil(maximum(abs.(error_stats.bias)), digits=1)
    rmse_diff_up_limit = ceil(maximum(abs.(Δrmse)), digits=1)
    colorranges = [
        (-rmse_diff_up_limit, rmse_diff_up_limit),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-bias_up_limit, bias_up_limit)
    ]

    n_metrics = length(x_labels)
    fig = Figure(size=fig_size)

    for (i, x_label) in enumerate(x_labels)
        yticks = i == 1 ?
                 string.(round.(y_coords[y_coords_sortperm], digits=2)) :
                 fill("", length(y_coords))
        ylabel = i == 1 ? "Latitude" : ""
        ax = Axis(
            fig[1, i],
            xlabel=x_label,
            ylabel=ylabel,
            xticks=([i], [""]),
            yticks=(1:length(y_coords), yticks),
        )
        hm = heatmap!(ax, fill(i, length(y_coords)), 1:length(y_coords), heat_data[:, i];
            colormap=colormaps[i], colorrange=colorranges[i])
        Colorbar(fig[end+1, 1:n_metrics], hm, vertical=false, label=x_label)
    end

    # TODO Add labels for south, center and north regions

    return fig
end
