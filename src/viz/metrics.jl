using CairoMakie
using GeoMakie

"""
    _plot_bootstrap_metric_scatter(estimate, ci_lo, ci_hi, block_eligible, median, median_lo,
                                    median_hi; title, ylabel, metric_label, fig_opts, axis_opts)::Figure

Shared scatter-plot core for per-reef bootstrap-CI metrics (as returned by
[`rmse_diff_stats`](@ref)/[`nse_stats`](@ref)): one dot per reef with a bootstrap CI
error bar (solid for block-bootstrap reefs with `n_years >= 5`, dashed for iid reefs with
`n_years == 4`), plus the aggregate median and its bootstrap CI (over
block-bootstrap-eligible reefs only) as a horizontal line/band. `title` is used verbatim
(callers build the metric-specific title text).

`y_margin`, if given, fixes the y-axis lower bound at `minimum(estimate) - y_margin` —
comfortably below the worst *point estimate*, but CI whiskers (which can be far wider,
e.g. a ratio-type metric like NSE blowing up for a near-zero-variance resample) are
truncated at that floor rather than stretching the axis to fit them. A truncated whisker
keeps its usual solid/dashed style (that already encodes block vs iid, so overloading it
to also mean "truncated" would be ambiguous) and gets a small downward-pointing marker
capping it at the floor instead, to signal "this CI continues further than shown." Leave
as `nothing` (the default) to auto-scale the lower bound to the data with no truncation —
appropriate for metrics whose CIs don't have this kind of unbounded-tail blowup. `y_high`,
if given, fixes the axis' upper bound (via `ylims!`) rather than auto-scaling — e.g. NSE
has a true upper bound of 1.
"""
function _plot_bootstrap_metric_scatter(
    estimate::Vector{Float64},
    ci_lo::Vector{Float64},
    ci_hi::Vector{Float64},
    block_eligible::BitVector,
    median_val::Float64,
    median_lo::Float64,
    median_hi::Float64;
    title::String,
    ylabel::String,
    metric_label::String,
    y_high::Union{Float64,Nothing}=nothing,
    y_margin::Union{Float64,Nothing}=nothing,
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    order = sortperm(estimate)
    estimate = estimate[order]
    ci_lo = ci_lo[order]
    ci_hi = ci_hi[order]
    block_eligible = block_eligible[order]

    y_low = isnothing(y_margin) ? nothing : minimum(estimate) - y_margin
    clipped_lo = isnothing(y_low) ? ci_lo : max.(ci_lo, y_low)
    truncated = isnothing(y_low) ? falses(length(ci_lo)) : ci_lo .< y_low

    axis_opts = merge(Dict{Symbol,Any}(:title => title, :ylabel => ylabel), axis_opts)

    size = get(fig_opts, :size, (600, 400))
    fig = Figure(; size=size)
    ax = Axis(
        fig[1, 1],
        title=axis_opts[:title],
        xlabel=get(axis_opts, :xlabel, "Index"),
        ylabel=axis_opts[:ylabel],
    )

    block_color = "#0072B2"
    iid_color = "#E69F00"
    median_color = "#009E73"

    x = collect(1:length(estimate))
    iid_mask = .!block_eligible
    trunc_block = truncated .& block_eligible
    trunc_iid = truncated .& iid_mask

    scatter!(ax, x[block_eligible], estimate[block_eligible], color=block_color)
    scatter!(ax, x[iid_mask], estimate[iid_mask], color=iid_color)

    rangebars!(
        ax, x[block_eligible], clipped_lo[block_eligible], ci_hi[block_eligible],
        color=block_color, linestyle=:solid
    )
    rangebars!(
        ax, x[iid_mask], clipped_lo[iid_mask], ci_hi[iid_mask],
        color=iid_color, linestyle=:dash
    )

    if !isnothing(y_low)
        scatter!(
            ax, x[trunc_block], fill(y_low, sum(trunc_block));
            marker=:dtriangle, color=block_color, markersize=10
        )
        scatter!(
            ax, x[trunc_iid], fill(y_low, sum(trunc_iid));
            marker=:dtriangle, color=iid_color, markersize=10
        )
    end

    hlines!(ax, [median_val], color=median_color, linestyle=:solid, linewidth=2)
    hlines!(ax, [median_lo, median_hi], color=(median_color, 0.6), linestyle=:dash)

    legend_els = [
        MarkerElement(color=block_color, marker=:circle),
        MarkerElement(color=iid_color, marker=:circle),
        LineElement(color=median_color, linestyle=:solid),
    ]
    legend_labels = [
        "$metric_label (block bootstrap, n≥5 yrs)",
        "$metric_label (iid bootstrap, n=4 yrs)",
        "Median [95% CI]",
    ]
    if any(truncated)
        push!(legend_els, MarkerElement(color=:gray, marker=:dtriangle))
        push!(legend_labels, "CI truncated at axis floor")
    end

    Legend(fig[1, 2], legend_els, legend_labels)

    ylims!(ax, y_low, y_high)

    return fig
end

"""
    plot_rmse_scatter(rmse_stats::NamedTuple; observation_type="", fig_opts=Dict(), axis_opts=Dict())::Figure

Scatter plot of per-reef `benchmark_RMSE - model_RMSE` (as returned by
[`rmse_diff_stats`](@ref)) — see [`_plot_bootstrap_metric_scatter`](@ref).
"""
function plot_rmse_scatter(
    rmse_stats::NamedTuple;
    observation_type::String="",
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    n_eligible = sum(rmse_stats.block_eligible)
    n_greater_than_zero = sum(rmse_stats.diff[rmse_stats.block_eligible] .> 0)
    pct_greater_than_zero = round(100 * n_greater_than_zero / n_eligible; digits=1)
    median_val = trunc(rmse_stats.median; digits=3)
    median_lo = trunc(rmse_stats.median_lo; digits=3)
    median_hi = trunc(rmse_stats.median_hi; digits=3)

    obs_title = isempty(observation_type) ? "" : "\n$(titlecase(observation_type)) data"
    title = "Benchmark - Model (RMSE)" * obs_title *
            "\n# > 0: $n_greater_than_zero ($pct_greater_than_zero%) | Median: $median_val [$median_lo, $median_hi]"

    return _plot_bootstrap_metric_scatter(
        rmse_stats.diff, rmse_stats.ci_lo, rmse_stats.ci_hi, rmse_stats.block_eligible,
        rmse_stats.median, rmse_stats.median_lo, rmse_stats.median_hi;
        title=title, ylabel="RMSE Difference", metric_label="RMSE Diff",
        fig_opts=fig_opts, axis_opts=axis_opts,
    )
end

"""
    plot_nse_scatter(nse_stats::NamedTuple; observation_type="", fig_opts=Dict(), axis_opts=Dict())::Figure

Scatter plot of per-reef Nash-Sutcliffe efficiency (as returned by
[`nse_stats`](@ref)) — see [`_plot_bootstrap_metric_scatter`](@ref). NSE normalizes
`benchmark_RMSE - model_RMSE` by each reef's own variance (`RMSE_benchmark`), so reefs
with little natural variability aren't penalized relative to more volatile ones purely
for having less absolute room to beat the benchmark.
"""
function plot_nse_scatter(
    nse_stats::NamedTuple;
    observation_type::String="",
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    n_eligible = sum(nse_stats.block_eligible)
    n_greater_than_zero = sum(nse_stats.nse[nse_stats.block_eligible] .> 0)
    pct_greater_than_zero = round(100 * n_greater_than_zero / n_eligible; digits=1)
    median_val = trunc(nse_stats.median; digits=3)
    median_lo = trunc(nse_stats.median_lo; digits=3)
    median_hi = trunc(nse_stats.median_hi; digits=3)

    obs_title = isempty(observation_type) ? "" : "\n$(titlecase(observation_type)) data"
    title = "Nash-Sutcliffe Efficiency" * obs_title *
            "\n# > 0: $n_greater_than_zero ($pct_greater_than_zero%) | Median: $median_val [$median_lo, $median_hi]"

    return _plot_bootstrap_metric_scatter(
        nse_stats.nse, nse_stats.ci_lo, nse_stats.ci_hi, nse_stats.block_eligible,
        nse_stats.median, nse_stats.median_lo, nse_stats.median_hi;
        title=title, ylabel="NSE", metric_label="NSE", y_high=1.0, y_margin=1.0,
        fig_opts=fig_opts, axis_opts=axis_opts,
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
