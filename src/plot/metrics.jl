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
    mean_rmse_diff = mean(_rmse_diff)
    std_rmse_diff = std(_rmse_diff)

    n_greater_than_zero = sum(_rmse_diff .> 0)
    success_rate = round((n_greater_than_zero / length(_rmse_diff)) * 100, digits=2)

    obs_title = isempty(observation_type) ? "" : "\n$(titlecase(observation_type)) data"
    axis_opts::Dict{Symbol,Any} = Dict(
        :title => "Benchmark - Model (RMSE)" * obs_title,
        :ylabel => "RMSE Difference"
    )

    opts::Dict{Symbol,Any} = Dict(:metric_label => "RMSE Diff",)

    return plot_metric_scatter(
        _rmse_diff, mean_rmse_diff;
        fig_opts=fig_opts, axis_opts=axis_opts, opts=opts
    )
end

function plot_pcc_scatter(
    pcc;
    observation_type::String="",
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    lower_bound, mean, upper_bound = average_cc(pcc)

    _yticks = sort(unique(vcat((-1:0.5:1), round(mean, digits=2))))
    obs_title = isempty(observation_type) ? "" : "\n$(titlecase(observation_type)) data"
    axis_opts::Dict{Symbol,Any} = Dict(
        :title => "Pearson Correlation Coefficient (PCC)" * obs_title,
        :ylabel => "PCC",
        :yticks => _yticks
    )

    opts::Dict{Symbol,Any} = Dict(:metric_label => "PCC",)

    return plot_metric_scatter(
        sort(pcc), mean;
        fig_opts=fig_opts, axis_opts=axis_opts, opts=opts
    )
end

function plot_metric_scatter(
    metric_data, agg_data;
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)
    size = get(fig_opts, :size, (600, 400))

    title = get(axis_opts, :title, "")
    xlabel = get(axis_opts, :xlabel, "Index")
    ylabel = get(axis_opts, :ylabel, "Metric value")
    yticks = get(axis_opts, :yticks, Makie.automatic)

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
    raw_data::Array{Float64,3};
    fig_opts::Dict=Dict(),
    observations::LocationDataStore=COMBINED_STORE,
)::Figure
    n_validation_obs = length(observations.ltmp_unique_ids)
    ltmp_loc_indexes = collect(1:n_validation_obs)

    error_stats = collect_error_stats.([raw_data], ltmp_loc_indexes; observations=observations)
    rmse_, benchmark_, cc_, maee_, bias_, scc_ = eachrow(hcat(map(collect, error_stats)...))
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

    observation_gpkg = domain_gpkg[[findfirst(domain_gpkg.UNIQUE_ID .== ltmp_id) for ltmp_id in observations.ltmp_unique_ids], :]

    try
        poly!.(ax, domain_gpkg.geom, color=:gray)
        poly!.(ax, observation_gpkg.geom, color=:red)
    catch
        poly!.(ax, domain_gpkg.geometry, color=:gray)
        poly!.(ax, observation_gpkg.geometry, color=:red)
    end

    max_val, min_val = extrema(rmse_diffs)
    up_limit = maximum(abs.((max_val, min_val)))
    lower_limit = -up_limit

    try
        scatter!(ax, observation_gpkg.X_COORD, observation_gpkg.Y_COORD;
            markersize=35, colorrange=(lower_limit, up_limit), color=rmse_diffs, colormap=:bam, alpha=0.8,
            strokewidth=1, strokecolor=(:gray, 0.1)
        )
    catch
        scatter!(ax, observation_gpkg.LON, observation_gpkg.LAT;
            markersize=35, colorrange=(lower_limit, up_limit), color=rmse_diffs, colormap=:bam, alpha=0.8,
            strokewidth=1, strokecolor=(:gray, 0.1)
        )
    end

    Colorbar(
        fig[1, 2];
        colorrange=(lower_limit, up_limit), colormap=:bam, label="Benchmark RMSE - Model RMSE"
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
    rmse_, benchmark_, cc_, maee_, bias_, scc_ = eachrow(hcat(map(collect, error_stats)...))

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
    )
    observations = VALIDATION_STORE
    domain_gpkg = observations.domain_gpkg
    domain_gpkg[[findfirst(domain_gpkg.UNIQUE_ID .== ltmp_id) for ltmp_id in observations.ltmp_unique_ids], :]
    observation_gpkg = domain_gpkg[[findfirst(domain_gpkg.UNIQUE_ID .== ltmp_id) for ltmp_id in observations.ltmp_unique_ids], :]# domain_gpkg[domain_gpkg.UNIQUE_ID.∈[observations.ltmp_unique_ids], :]
    try
        poly!.(ax, domain_gpkg.geom, color=:gray)
    catch
        poly!.(ax, domain_gpkg.geometry, color=:gray)
    end

    max_val, min_val = extrema(cc_)
    up_limit = maximum(abs.((max_val, min_val)))
    lower_limit = -up_limit

    try
        scatter!(ax, observation_gpkg.X_COORD, observation_gpkg.Y_COORD, markersize=35,
            color=cc_, colorrange=(lower_limit, up_limit), colormap=:bam, alpha=0.6, strokewidth=1, strokecolor=(:gray, 0.1))
    catch
        scatter!(ax, observation_gpkg.LON, observation_gpkg.LAT, markersize=35,
            color=cc_, colorrange=(lower_limit, up_limit), colormap=:bam, alpha=0.6, strokewidth=1, strokecolor=(:gray, 0.1))
    end
    Colorbar(fig[1, 2], colorrange=(lower_limit, up_limit), colormap=:bam,
        label="Pearson Correlation Coefficient")

    return fig
end

function plot_metrics_heatmap(rs_raw; fig_size=(500, 700), observations=COMBINED_STORE)::Figure
    n_validation_obs = length(observations.ltmp_unique_ids)
    ltmp_loc_indexes = collect(1:n_validation_obs)
    error_stats = collect_error_stats.(
        Ref(rs_raw), ltmp_loc_indexes; observations=observations
    )
    rmse_, benchmark_, cc_, maee_, bias_, scc_ = eachrow(hcat(map(collect, error_stats)...))

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
    heat_data = hcat(benchmark_ .- rmse_, cc_, scc_, bias_)[y_coords_sortperm, :]
    x_labels = ["Benchmark - Model (RMSE)", "PCC", "SCC", "Bias"]
    colormaps = [:Spectral, :Spectral, :Spectral, :Spectral]

    bias_up_limit = ceil(maximum(abs.(bias_)), digits=1)
    rmse_diff_up_limit = ceil(maximum(abs.(benchmark_ .- rmse_)), digits=1)
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

    return fig
end
