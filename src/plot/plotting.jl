# TODO Rename/organize this file

"""
    plot_class_size_props(init_cover, class_idx)::Figure
"""
function plot_class_size_props(
    init_cover,
    class_idx
)::Figure
    class_state = init_cover[(class_idx-1)*11+1:class_idx*11]
    taxa_names = ADRIA.functional_group_names()
    taxa_prop = class_state[7:11]
    f = Figure(; size=(1200, 900))
    ax = Axis(f[1, 1], xlabel="taxonomy", xticks=(1:5, String.(taxa_names)), ylabel="taxa size lambda", title="class: $(class_idx)")
    barplot!(
        ax,
        1:5,
        taxa_prop
    )
    return f
end

"""
    plot_class_properties(init_cover, class_idx)::Figure
"""
function plot_class_properties(
    init_cover,
    class_idx
)::Figure
    class_state = init_cover[(class_idx-1)*11+1:class_idx*11]
    taxa_names = ADRIA.functional_group_names()
    taxa_prop = class_state[2:6] ./ sum(class_state[2:6])
    f = Figure(; size=(1200, 900))
    ax = Axis(f[1, 1], xlabel="taxonomy", xticks=(1:5, String.(taxa_names)), ylabel="taxa proportions", title="class: $(class_idx)")
    barplot!(
        ax,
        1:5,
        taxa_prop
    )
    return f
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

    fig_size = get(fig_opts, :size, (500, 700))
    fig = Figure(size=fig_size)
    ax = Axis(fig[1, 1], title="Benchmark RMSE - Model RMSE", titlesize=18, xlabelsize=15,
        ylabelsize=15, xlabel="Longitude", ylabel="Latitude")
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

function plot_pearson_coeff_map(
    raw_data::Array{Float64,3};
    fig_opts::Dict=Dict(),
    observations::LocationDataStore=COMBINED_STORE,
)::Figure
    raw_data = rs_raw.raw
    n_validation_obs = length(observations.ltmp_unique_ids)
    ltmp_loc_indexes = collect(1:n_validation_obs)

    error_stats = collect_error_stats.([raw_data], ltmp_loc_indexes; observations=observations)
    rmse_, benchmark_, cc_, maee_, bias_ = eachrow(hcat(map(collect, error_stats)...))

    fig_size = get(fig_opts, :size, (500, 700))
    fig = Figure(size=fig_size)
    ax = Axis(fig[1, 1], title="Pearson's Correlation Coefficient", titlesize=18, xlabelsize=15,
        ylabelsize=15, xlabel="Longitude", ylabel="Latitude")
    domain_gpkg = observations.domain_gpkg
    observation_gpkg = domain_gpkg[domain_gpkg.UNIQUE_ID.∈[observations.ltmp_unique_ids], :]
    poly!.(ax, domain_gpkg.geom, color=:gray)
    poly!.(ax, observation_gpkg.geom, color=:red)
    scatter!(ax, observation_gpkg.X_COORD, observation_gpkg.Y_COORD, markersize=35,
        color=cc_, colormap=:bam, alpha=0.8, strokewidth=1, strokecolor=(:gray, 0.1))
    max_val, min_val = extrema(cc_)
    up_limit = maximum(abs.((max_val, min_val)))
    lower_limit = -up_limit
    Colorbar(fig[1, 2], limits=(lower_limit, up_limit), colormap=:bam,
        label="Model RMSE - Benchmark RMSE")

    return fig
end

function plot_ltmp_disturbances!(
    fig::Union{Figure,GridLayout}, ax_row::Int64, loc_disturbances; axis_opts::Dict=Dict()
)::Nothing
    limits = get(axis_opts, :limits, nothing)
    xlimits = !isnothing(limits) ? limits[1] : nothing
    axis_opts[:limits] = (xlimits, (0, 1))
    #@info(axis_opts[:limits])
    axis_opts = Dict(
        axis_opts...,
        :yticks => [0, 1],
        :yticklabelsize => 0,
        :ylabel => "Disturbances",
    )

    ax_ltmp = Axis(fig[ax_row, 1]; axis_opts...)


    disturbance_vlines, disturbance_types = plot_ltmp_disturbances!(ax_ltmp, loc_disturbances)

    return nothing
end
function plot_ltmp_disturbances!(ax::Axis, loc_disturbances)
    disturbances_types = collect(loc_disturbances.disturbances)
    n_disturbances_types = length(disturbances_types)
    disturbances_years = collect(loc_disturbances.timesteps)
    plotted_vlines = []
    plotted_types = []
    for (disturbance_type_idx, disturbance_type) in enumerate(disturbances_types)
        disturbance_years_mask = loc_disturbances[disturbances=At(disturbance_type)] .== 1
        if any(disturbance_years_mask)
            vline = vlines!(
                ax,
                disturbances_years[disturbance_years_mask],
                color=disturbance_type_idx,
                linewidth=3,
                colormap=COLORMAPS.ltmp_disturbances,
                colorrange=(1, n_disturbances_types)
            )
            push!(plotted_vlines, vline)
            push!(plotted_types, disturbance_type)
        end
    end

    return plotted_vlines, plotted_types
end

"""
    taxa_cover_proportions(raw_data)::Figure

Plot the proportion of coral cover composed of each functional group as a line graph given
the raw modelled output cover matrix as input.
"""
function taxa_cover_proportions(raw_data)::Figure
    cover = reshape(raw_data, (15, 7, 5, 3806))
    cover = dropdims(sum(cover, dims=4), dims=4)
    cover ./= sum(cover, dims=(2, 3))
    cover = dropdims(sum(cover, dims=2), dims=2)
    cover = permutedims(cover, (2, 1))
    xs = 2008:2022

    f = Figure(; size=(1200, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Proportion",
        title="Functional Group Cover Proportions",
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, cover, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    # Legend(f[1, 2], ax, framevisible=false)
    return f
end

"""
    taxa_population_proportions(raw_data)::Figure

Plot the proportion of coral population composed of each functional group. Population
proportions are calculated using the average coral diameter of each size class.
"""
function taxa_population_proportions(raw_data)::Figure
    sc_mean_area = reshape(
        permutedims(ADRIA.colony_areas()[2], (2, 1)),
        (1, 35, 1)
    )
    population = raw_data ./ sc_mean_area
    population = reshape(population, (15, 7, 5, 3806))
    population = dropdims(sum(population, dims=4), dims=4)
    population ./= sum(population, dims=(2, 3))
    population = dropdims(sum(population, dims=2), dims=2)
    population = permutedims(population, (2, 1))
    xs = 2008:2022

    f = Figure(; size=(1200, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title="Functional Group Population Proportions",
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, population, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    #Legend(f[1, 2], ax, framevisible=false)
    return f
end

"""
    temporal_size_class_proportions(raw_data)

Calculate the percentage of coral population occupied by each size class split by functional
group.
"""
function temporal_size_class_proportions(raw_data)::Figure
    sc_mean_area = reshape(
        permutedims(ADRIA.colony_areas()[2], (2, 1)),
        (1, 35, 1)
    )
    n_locs = size(raw_data, 3)
    population = raw_data ./ sc_mean_area
    population = reshape(population, (15, 7, 5, n_locs))
    population = dropdims(sum(population, dims=4), dims=4)
    population ./= sum(population, dims=2)
    population = permutedims(population, (3, 2, 1))

    fg_names = ADRIA.functional_group_names()
    xs = 2008:2022
    col = :oslo10

    f = Figure(; size=(1600, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[1]),
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, population[1, :, :], color=col, labels="Size Class: " .* string.(1:7))
    Axis(
        f[1, 2];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[2]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[2, :, :], color=col)
    Axis(
        f[2, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[3]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[3, :, :], color=col)
    Axis(
        f[2, 2];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[4]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[4, :, :], color=col)
    Axis(
        f[1, 3];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[5]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[5, :, :], color=col)
    Legend(f[2, 3], ax, framevisible=false)
    colsize!(f.layout, 3, Relative(1 / 3))
    resize_to_layout!(f)
    return f
end

function plot_ltmp(ltmp_n, ltmp_c, ltmp_s)::Nothing
    f = Figure(; size=(900, 1600))
    Axis(
        f[1, 1],
        title="North GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )

    band!(ltmp_n.Year, ltmp_n.lower, ltmp_n.upper, color=(:black, 0.4))
    lines!(ltmp_n.Year, ltmp_n.response, color=:black)
    Axis(
        f[2, 1],
        title="Central GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )


    band!(ltmp_c.Year, ltmp_c.lower, ltmp_c.upper, color=(:black, 0.4))
    lines!(ltmp_c.Year, ltmp_c.response, color=:black)
    Axis(
        f[3, 1],
        title="South GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )

    band!(ltmp_s.Year, ltmp_s.lower, ltmp_s.upper, color=(:black, 0.4))
    lines!(ltmp_s.Year, ltmp_s.response, color=:black)
    save("LTMP_data.png", f)
    return nothing
end

"""
    plot_all_regions(dom, model_results; region_masks=[NORTH_MASK, CENTRAL_MASK, SOUTH_MASK], ref_years=START_YEAR:END_YEAR)::Figure
    plot_all_regions(north_results, north_obs, central_results, central_obs, south_results, south_obs)::Figure

Plot the modelled results against observations/modelled observation results.
"""
function plot_all_regions(
    dom,
    model_results;
    region_masks=[NORTH_MASK, CENTRAL_MASK, SOUTH_MASK],
    ref_years=START_YEAR:END_YEAR
)::Figure
    s_rac = (dropdims(sum(model_results.raw, dims=2), dims=2) .* site_k_area(dom)') ./ loc_area(dom)'

    ref_years = START_YEAR:END_YEAR

    north_res = ADRIA.DataCube(s_rac[:, NORTH_MASK, :]; timesteps=ref_years, sites=1:count(NORTH_MASK), scenarios=1:1)
    central_res = ADRIA.DataCube(s_rac[:, CENTRAL_MASK, :]; timesteps=ref_years, sites=1:count(CENTRAL_MASK), scenarios=1:1)
    south_res = ADRIA.DataCube(s_rac[:, SOUTH_MASK, :]; timesteps=ref_years, sites=1:count(SOUTH_MASK), scenarios=1:1)

    f = plot_all_regions(
        north_res, ltmp_north, central_res, ltmp_central, south_res, ltmp_south
    )

    return f
end
function plot_all_regions(
    north_results, north_obs, central_results, central_obs, south_results, south_obs
)::Figure

    f = Figure(; size=(1600, 1600))
    north_ax = plot_region(
        f,
        1,
        1,
        "North GBR",
        north_obs,
        north_results
    )
    central_ax = plot_region(
        f,
        1,
        2,
        "Central GBR",
        central_obs,
        central_results
    )
    south_ax = plot_region(
        f,
        2,
        1,
        "South GBR",
        south_obs,
        south_results;
        showlegend=true,
        legend_row=2,
        legend_col=2
    )
    linkyaxes!(north_ax, central_ax, south_ax)

    resize_to_layout!(f)
    return f
end

"""
    plot_region(f::Figure, row::Int64, col::Int64, title::String, obs::DataFrame, sim::YAXArray; showlegend::Bool = false, legend_row::Int64 = 2, legend_col::Int64 = 2)::Axis

"""
function plot_region(
    f::Figure,
    row::Int64,
    col::Int64,
    title::String,
    obs::DataFrame,
    sim::YAXArray;
    showlegend::Bool=false,
    legend_row::Int64=2,
    legend_col::Int64=2
)::Axis

    ax = Axis(
        f[row, col],
        title=title,
        xlabel="year",
        ylabel="Relative Absolute Cover",
        width=400,
        height=400
    )

    mean_agg = dropdims(mean(sim, dims=:sites), dims=:sites)
    confints = ADRIA.analysis.series_confint(mean_agg.data)

    xs::Vector{Float64} = collect(sim.timesteps)
    ADRIA_series = series!(xs, sim[:, :, 1].data[:, :]'; solid_color=(:red, 0.01), linewidth=1, labels=nothing)
    ADRIA_line = lines!(xs, confints[:, 2], color=:red, linewidth=5)

    xs = obs.Year
    lower::Vector{Float64} = obs.lower
    upper::Vector{Float64} = obs.upper
    response::Vector{Float64} = obs.response

    LTMP_band = band!(xs, lower, upper, color=(:black, 0.4))
    LTMP_line = lines!(xs, response, color=:black)

    if showlegend
        Legend(
            f[legend_row, legend_col],
            [[LTMP_line, LTMP_band], [ADRIA_line, ADRIA_series]],
            ["LTMP", "CoralBlox"]
        )
    end

    return ax
end

function plot_residual(
    f::Figure,
    ax_row::Int64,
    ax_col::Int64,
    ltmp_data::DataFrame,
    ADRIA_data::YAXArray,
    title::String
)::Nothing

    Axis(
        f[ax_row, ax_col],
        title=title,
        xlabel="year",
        ylabel="Residuals"
    )

    ltmp_ind::Int64 = findfirst(x -> x >= START_YEAR, ltmp_data.Year)
    ltmp_end::Int64 = findfirst(x -> x >= END_YEAR, ltmp_data.Year) - 1

    ltmp_xs = ltmp_data.Year[ltmp_ind:ltmp_end]
    ltmp_resp = ltmp_data.response[ltmp_ind:ltmp_end]

    adria_mean = dropdims(mean(
            ADRIA_data[timesteps=At(ltmp_xs)], dims=(:sites, :scenarios)
        ), dims=(:sites, :scenarios))

    residuals::Vector{Float64} = ltmp_resp .- adria_mean
    scatter!(ltmp_xs, residuals, color=:black)

    return nothing
end

function plot_residuals(
    filename::String; fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:size => (1600, 1600))
)::Figure

    f = Figure(; fig_opts...)

    plot_residual(
        f,
        1,
        1,
        ltmp_north,
        north_res,
        "North GBR"
    )

    plot_residual(
        f,
        2,
        1,
        ltmp_central,
        central_res,
        "Central GBR"
    )

    plot_residual(
        f,
        3,
        1,
        ltmp_south,
        south_res,
        "South GBR"
    )

    save(filename, f)

    return f
end

function create_compare_plot(
    filename::String; fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:size => (1600, 1600))
)::Figure

    f = Figure(; fig_opts...)

    plot_region(
        f,
        1,
        1,
        "North GBR",
        ltmp_north,
        north_res
    )
    plot_region(
        f,
        1,
        2,
        "Central GBR",
        ltmp_central,
        central_res
    )
    plot_region(
        f,
        2,
        1,
        "South GBR",
        ltmp_south,
        south_res;
        showlegend=true,
        legend_row=2,
        legend_col=2
    )
    resize_to_layout!(f)

    save(filename, f)
    return f
end

"""
    plot_coral_param(loc::String, param_name::String, category, group, data)::Figure

Plot coral parameters in a  bar graph grouped by size class.

# Example

```julia
# Each parameter size class. [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, ...]
category = [floor((i - 1) / 5) + 1 for i in 1:35]
# Each parameter functional group, [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, ...]
group = [(i - 1) % 5 + 1 for i in 1:35]
# Flattened coral background mortality rate.
flt_mb_rate = corals.mb_rate
# Plot parameters
f = plot_coral_param(<Loc Name>, "Background Mortality", category, group, flt_mb_rate)
```

# Arguments
- `loc` : Name of location. Used in title.
- `param_name` : Name of parameter being plotted. Used in title.
- `category` : A vector of ints describing the size class each parameter belongs to.
- `group` : A vector of ints describing the functional groups each parameter belongs to.
- `data` : A vector of parameters to be plotted.
"""
function plot_coral_param(
    loc::String,
    param_name::String,
    category,
    group,
    data
)::Figure

    fig = Figure(; size=(1300, 900))
    ax = Axis(
        fig[1, 1],
        xticks=1:length(category)/5,
        title="$(loc): $(param_name)"
    )

    barplot!(
        ax,
        category,
        data,
        dodge=group,
        color=Makie.wong_colors()[group]
    )

    labels = String.(ADRIA.functional_group_names())
    elements = [PolyElement(polycolor=Makie.wong_colors()[i]) for i in 1:length(labels)]
    title = "Functional Groups"

    Legend(fig[1, 2], elements, labels, title)

    return fig
end

"""
    plot_taxa_props(loc::String, cover)::Figure

Plot the modelled coral composition for a given location.

# Example
```julia
# ADRIA single run results
rs = ADRIA.run_model(...)
# Index of location of interest
location_idx = ...

f = plot_taxa_props(<Location Name or ID>, rs.raw[:, :, location_idx])
```
"""
function plot_taxa_props(
    loc::String,
    cover;
    opts::Dict=Dict()
)::Figure
    f = Figure(; size=(1300, 900))
    ax_opts = Dict(
        :xlabel => "year",
        :ylabel => "Cover Proportion",
        :title => "$(loc): Functional Group Cover Proportions",
        :limits => (nothing, nothing, 0, 1)
    )
    plot_taxa_props!(f, 1, cover; ax_opts=ax_opts, opts=opts)
    return f
end
function plot_taxa_props!(
    fig::Union{Figure,GridLayout},
    ax_row::Int64,
    cover;
    ax_opts::Dict=Dict(),
    opts::Dict=Dict()
)::Nothing
    ax_opts = Dict(ax_opts..., :xlabel => "Year", :ylabel => "Proportional Cover",)
    ax_taxa = Axis(fig[ax_row, 1]; ax_opts...)
    taxa_plot = plot_taxa_props!(ax_taxa, cover; opts=opts)
    # Legend(fig[ax_row, 2], taxa_plot)
    nothing
end
function plot_taxa_props!(
    ax::Axis,
    cover;
    opts::Dict=Dict()
)::Axis
    color = COLORMAPS.taxa
    cover = reshape(cover, (15, 7, 5))
    cover = dropdims(sum(cover, dims=2), dims=2) ./ dropdims(sum(cover, dims=(2, 3)), dims=2)
    cover = permutedims(cover, (2, 1))
    xs = 2008:2022
    series!(
        xs,
        cover,
        color=color,
        labels=String.(ADRIA.functional_group_names()),
        linewidth=3,
    )
    return ax
end

"""
    plot_bleaching_mortality(loc::String, loc_cover, loc_bleaching)::Figure

Plot the percentage of cover lost to bleaching over time.

# Example
```julia
# ADRIA single run results
rs = ADRIA.run_model(...)
# Index of location of interest
location_idx = ...

f = plot_bleaching_mortality(
    <Location Name or ID>, # Only used in title.
    rs.raw[:, :, location_idx],
    rs.bleaching_mortality[:, :, location_idx]
)
```
"""
function plot_bleaching_mortality(
    loc::String,
    loc_cover,
    loc_bleaching
)::Figure
    prop_sc_fg_cover = permutedims(reshape(loc_cover, (15, 7, 5)), (1, 3, 2))
    prop_sc_fg_cover = prop_sc_fg_cover ./ sum(prop_sc_fg_cover, dims=(2, 3))
    perc_loss = dropdims(sum(prop_sc_fg_cover .* loc_bleaching, dims=(2, 3)), dims=(2, 3))

    xs = START_YEAR:END_YEAR
    f = Figure(; size=(1300, 900))
    Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Loss",
        title="$(loc): Proportion Cover Loss from Bleaching"
    )
    lines!(xs, perc_loss)

    return f
end

"""
    plot_cyclone_mortality(loc::String, loc_cover, loc_cyclone)::Figure

Plot the percentage of cover lost to cyclones over time.

# Example
```julia
# ADRIA single run results
rs = ADRIA.run_model(...)
# Index of location of interest
location_idx = ...
# Get cyclone scenario used
cyc_scen = scens[1, :cyclone_mortality_scenario]

f = plot_cyclone_mortality(
    <Location Name or ID>, # Only used in title.
    rs.raw[:, :, location_idx],
    dom.cyclone_mortality_scens[:, location_idx, :, cyc_scen].data
)
```
"""
function plot_cyclone_mortality(
    loc::String,
    loc_cover,
    loc_cyclone
)::Figure
    prop_sc_fg_cover = permutedims(reshape(loc_cover, (15, 7, 5)), (1, 3, 2))
    prop_sc_fg_cover = dropdims(sum(prop_sc_fg_cover ./ sum(prop_sc_fg_cover, dims=(2, 3)), dims=3), dims=3)
    perc_loss = dropdims(sum(prop_sc_fg_cover .* loc_cyclone, dims=2), dims=2)

    xs = 2008:2022
    f = Figure(; size=(1300, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Loss",
        title="$(loc): Proportion Cover Loss from Cyclones"
    )
    sr = lines!(xs, perc_loss)
    return f
end

"""
    plot_target_proportions(loc::String, target_taxa)::Figure

Plot the target coral composition over time.

# Example
```julia
f = plot_target_proportions(<Location Name or ID>, rm_ltmp_taxa[:, :, limited_loc_pos])
```
"""
function plot_target_proportions(loc::String, target_taxa)::Figure
    xs = 2008:2022
    normalised_comp = target_taxa ./ sum(target_taxa, dims=2)
    non_missing_mask = (!).(ismissing.(target_taxa[:, 1]))

    f = Figure(; size=(1300, 900))

    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Loss",
        xticks=xs,
        title="$(loc): Target Coral Composition",
        limits=(nothing, nothing, 0, 1)
    )

    xs_m = xs[non_missing_mask]
    dt = permutedims(Float64.(normalised_comp[non_missing_mask, :]), (2, 1))

    sr = series!(xs_m, dt, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    Legend(f[1, 2], ax, framevisible=false)
    return f
end

function _mortality_to_cyc_category(cyc_mortality_scens)
    for (idx_mort, mortality) in enumerate([0.0, 0.0121482, 0.0155069, 0.0197484, 0.0246357, 0.0302982])
        cyc_mortality_scens[cyc_mortality_scens.==mortality] .= idx_mort - 1
    end
    return cyc_mortality_scens
end
