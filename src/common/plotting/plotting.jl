include("plot_helpers.jl")

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

function plot_modelled_v_ltmp!(
    fig::Figure,
    ax_row::Int64,
    raw_data,
    ltmp_loc_idx;
    observations::LocationDataStore=COMBINED_STORE,
    dom=dom,
    axis_opts::Dict=Dict()
)::Nothing
    ax = Axis(fig[ax_row, 1]; axis_opts...)
    obs, sim = plot_modelled_v_ltmp!(ax, raw_data, ltmp_loc_idx; observations=observations)
    Legend(f[ax_row, 2], [obs, sim], ["LTMP", "CoralBlox"], "Coral Cover")
    return nothing
end
function plot_modelled_v_ltmp!(
    ax::Axis,
    raw_data,
    ltmp_loc_idx;
    observations::LocationDataStore=COMBINED_STORE,
    dom=dom
)
    domain_loc_idx::Int64 = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    loc_k_area = ADRIA.site_k_area(dom)[domain_loc_idx]
    loc_area = ADRIA.loc_area(dom)[domain_loc_idx]

    loc_cover = dropdims(
        sum(raw_data[:, :, domain_loc_idx], dims=2),
        dims=2) .* loc_k_area ./ loc_area

    obs_loc_data = observations.ltmp_coral_cover[ltmp_loc_idx, :]
    not_missing_obs = (!).(ismissing.(obs_loc_data))
    obs_tf = (START_YEAR:END_YEAR)[not_missing_obs]

    obs = scatter!(ax, obs_tf, obs_loc_data[not_missing_obs], color=(:red, 0.5), markersize=20)
    sim = lines!(ax, 2008:2022, loc_cover, color=:red)
    return obs, sim
end

"""
    location_comparison( raw_data,  ltmp_loc_idx,  save_dir;  obs_idxs=ALL_LTMP_IDXS,  obs_loc_labels=ltmp_reef_data.RME_UNIQUE_ID)::Figure
Plot LTMP Manta Tow Coral Cover against the modelled cover output given the LTMP location index.
"""
function location_comparison(
    dom,
    raw_data,
    ltmp_loc_idx,
    cyc_scens,
    dhw_scens,
    disturbances;
    observations::LocationDataStore=COMBINED_STORE
)::Figure
    # Extra information specific to the given location
    reef_id = observations.ltmp_unique_ids[ltmp_loc_idx]

    f = Figure(; size=(1000, 1000))

    title_text = _location_err_title(
        raw_data,
        ltmp_loc_idx;
        observations=observations
    )

    start_year, end_year = extrema(cyc_scens.timesteps)
    xlimits = (start_year - 0.5, end_year + 0.5)
    xticks = start_year:end_year

    cover_limits = (xlimits, nothing)
    base_axis_opts::Dict = Dict(
        :xlabelsize => 20,
        :ylabelsize => 20,
        :xticks => xticks,
        :xticklabelrotation => pi / 6,
        :limits => cover_limits,
    )

    # Plot LTMP and CoralBlox covers
    # cover_limits = (xlimits, (0, maximum(raw_data) * 1.1))
    modelled_v_ltmp_row = 1
    modelled_ltmp_axis_opts = Dict(
        base_axis_opts..., :title => title_text, :ylabel => "Relative total cover",
    )
    plot_modelled_v_ltmp!(
        f, modelled_v_ltmp_row, raw_data, ltmp_loc_idx;
        observations=observations, axis_opts=modelled_ltmp_axis_opts
    )

    # Plot CoralBlox Observed DHW and Cyclone Categories
    ax_opts_coralblox_disturbances = Dict(base_axis_opts..., :ylabel => "DHW", :height => 100)
    ax_row_coralblox_disturbances = 2
    loc_dhw_scens = dhw_scens[locs=At(reef_id)]
    loc_cyc_scens = cyc_scens[locations=At(reef_id)]
    plot_coralblox_disturbances!(
        f, ax_row_coralblox_disturbances, loc_dhw_scens, loc_cyc_scens;
        axis_opts=ax_opts_coralblox_disturbances,
    )

    if reef_id ∈ disturbances.locations
        loc_ltmp_disturbances = disturbances[locations=At(reef_id)]
        ax_row_ltmp_disturbances = 3
        plot_ltmp_disturbances!(
            f, ax_row_ltmp_disturbances, loc_ltmp_disturbances;
            axis_opts=ax_opts_coralblox_disturbances
        )
        ax_row_taxa = 4
    else
        @warn("Reef $reef_id / $ltmp_loc_idx not found in LTMP disturbances DataFrame")
        ax_row_taxa = 3
    end

    cb_loc_id = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    cover = rs_raw.raw[:, :, cb_loc_id]
    plot_taxa_props!(f, ax_row_taxa, cover; ax_opts=base_axis_opts)
    return f
end

function plot_coralblox_disturbances!(
    fig::Figure, ax_row::Int64, loc_dhw_scens, loc_cyclone_scens; axis_opts::Dict=Dict()
)::Nothing
    ax_dhw = Axis(fig[ax_row, 1]; axis_opts...)
    dhw_plot = lines!(
        ax_dhw,
        parse.(Int64, loc_dhw_scens.timesteps),
        collect(loc_dhw_scens.data),
        color=:orange,
        linestyle=:dot,
        linewidth=4
    )

    ax_cyclone = Axis(fig[ax_row, 1]; axis_opts...)
    hidespines!(ax_cyclone)
    hidedecorations!(ax_cyclone)
    plotted_vlines, plotted_scens = plot_cyclone_scens!(ax_cyclone, loc_cyclone_scens)

    Legend(
        fig[ax_row, 2],
        [dhw_plot, plotted_vlines...], ["Observed DHW",
            plotted_scens...],
        "CoralBlox Disturbances"
    )
    return nothing
end

function plot_cyclone_scens!(ax::Axis, cyc_scens)
    unique_scens = sort(unique(cyc_scens)[unique(cyc_scens).>0])
    plotted_cyc_vlines = []
    plotted_cyc_scens = []
    for scen in unique_scens
        target_years = collect(cyc_scens[cyc_scens.==scen].timesteps)
        vline = vlines!(
            ax,
            target_years,
            color=Int64(scen + 1),
            colormap=(:grayC10),
            colorrange=(5, 1),
            linestyle=:dash,
            linewidth=3,
            ymin=0
        )
        push!(plotted_cyc_vlines, vline)
        push!(plotted_cyc_scens, "Cyc. Category: $(Int64(scen))")
    end
    return plotted_cyc_vlines, plotted_cyc_scens
end

function plot_ltmp_disturbances!(
    fig::Figure, ax_row::Int64, loc_disturbances; axis_opts::Dict=Dict()
)::Nothing
    limits = get(axis_opts, :limits, nothing)
    xlimits = !isnothing(limits) ? limits[1] : nothing
    axis_opts[:limits] = (xlimits, (0, 1))
    #@info(axis_opts[:limits])
    axis_opts = Dict(axis_opts..., :yticks => [0, 1], :yticklabelsize => 0, :height => 100)

    ax_ltmp = Axis(fig[ax_row, 1]; axis_opts...)


    disturbance_vlines, disturbance_types = plot_ltmp_disturbances!(ax_ltmp, loc_disturbances)
    Legend(fig[ax_row, 2], [disturbance_vlines...], [disturbance_types...], "LTMP Disturbances")
    return nothing
end
function plot_ltmp_disturbances!(ax::Axis, loc_disturbances; colormap=:seaborn_bright6)
    disturbances_types = collect(loc_disturbances.disturbances)
    n_disturbances_types = length(disturbances_types)
    colorrange = (1, n_disturbances_types)

    disturbances_years = collect(loc_disturbances.timesteps)
    plotted_vlines = []
    plotted_types = []
    _colors = ["#1f77b4", "#2ca02c", "#9467bd", "#17becf", "#bcbd22", "#e377c2"]
    for (disturbance_type_idx, disturbance_type) in enumerate(disturbances_types)
        disturbance_years_mask = loc_disturbances[disturbances=At(disturbance_type)] .== 1
        if any(disturbance_years_mask)
            vline = vlines!(
                ax,
                disturbances_years[disturbance_years_mask],
                color=disturbance_type_idx,
                linewidth=3,
                colormap=colormap,
                colorrange=colorrange
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
    Legend(f[1, 2], ax, framevisible=false)
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
    Legend(f[1, 2], ax, framevisible=false)
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
    cover
)::Figure
    f = Figure(; size=(1300, 900))
    ax_opts = Dict(
        :xlabel => "year",
        :ylabel => "Cover Proportion",
        :title => "$(loc): Functional Group Cover Proportions",
        :limits => (nothing, nothing, 0, 1)
    )
    plot_taxa_props!(f, 1, cover; ax_opts=ax_opts)
    return f
end
function plot_taxa_props!(
    fig::Figure,
    ax_row::Int64,
    cover;
    ax_opts::Dict=Dict()
)::Nothing
    ax_opts = Dict(ax_opts..., :xlabel => "Year", :ylabel => "Prop. Coral Composition",)
    ax_taxa = Axis(fig[ax_row, 1]; ax_opts...)
    taxa_plot = plot_taxa_props!(ax_taxa, cover)
    Legend(fig[ax_row, 2], taxa_plot)
    nothing
end
function plot_taxa_props!(
    ax::Axis,
    cover
)::Axis
    cover = reshape(cover, (15, 7, 5))
    cover = dropdims(sum(cover, dims=2), dims=2) ./ dropdims(sum(cover, dims=(2, 3)), dims=2)
    cover = permutedims(cover, (2, 1))
    xs = 2008:2022
    series!(xs, cover, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
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
