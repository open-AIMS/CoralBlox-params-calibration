"""
# Example
```
management_areas = unique(dom.loc_data.management_area_short)
management_area_masks = [dom.loc_data.management_area_short .== area for area in management_areas]

management_area_stats = region_stats(Symbol.(management_area_masks), management_areas)

plot_regional_comparison(management_area_stats)
```
"""
function plot_regional_comparison!(
    ax::Axis,
    model_stats::Matrix{Float64},
    historic_stats::Matrix{Float64},
    opts::Dict{Symbol,Any};
    xdata=START_YEAR:END_YEAR
)::Nothing
    band!(
        ax, xdata, model_stats[1, :], model_stats[3, :];
        color=(get(opts, :model_band_color, COLORS[:model_band_color]))
    )
    band!(
        ax, xdata, historic_stats[1, :], historic_stats[3, :];
        color=get(opts, :historic_band_color, COLORS[:historic_band_color])
    )

    scatterlines!(
        ax, xdata, collect(model_stats[2, :]);
        color=get(opts, :model_line_color, COLORS[:model_line_color]),
        markersize=get(opts, :markersize, 20)
    )
    scatterlines!(
        ax, xdata, historic_stats[2, :];
        color=get(opts, :historic_line_color, COLORS[:historic_line_color]),
        markersize=get(opts, :markersize, 20)
    )

    return nothing
end
function plot_regional_comparison(
    regional_stats::NamedTuple;
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    opts::Dict{Symbol,Any} = opts_regional_comparison(:blue, :red)
    f = Figure(; fig_opts...)
    plot_regional_comparison!(f, regional_stats; axis_opts=axis_opts, opts=opts)
    legend_regional_comparison!(f[end+1, :]; opts=opts)
    return f
end
function plot_regional_comparison!(
    grid::Union{GridPosition,GridLayout,Figure},
    regional_stats::NamedTuple;
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Nothing
    r_keys = keys(regional_stats)
    n_plots = length(r_keys)

    n_locations::Int64 = 0
    _rmse::Float64 = 0.0
    _mae::Float64 = 0.0
    for (idx, key) in enumerate(r_keys)
        n_locations = regional_stats[key].n_locations
        _rmse = trunc(regional_stats[key].rmse, digits=2)
        _mae = trunc(regional_stats[key].mae, digits=2)

        ax = Axis(
            grid[fig_coord(idx, n_plots)...],
            title="$(string(key)) ($n_locations Reefs)\nRMSE: $_rmse | MAE: $_mae",
            xlabel="Year",
            ylabel="Relative cover",
            # limits=axis_limits,
            xticks=START_YEAR:2:END_YEAR,
            xticklabelrotation=π / 4
        )

        plot_regional_comparison!(
            ax,
            regional_stats[key].model_stats,
            replace(regional_stats[key].historic_stats, missing => NaN),
            opts,
        )
    end

    Label(grid[0, :], get(opts, :title, ""); fontsize=18)

    return nothing
end
function plot_regional_comparison(
    regional_stats_validation::NamedTuple,
    regional_stats_calibration::NamedTuple;
    fig_size=(800, 800), axis_limits=(nothing, (0.0, 1.0)), fig_title=""
)::Figure
    opts_validation::Dict{Symbol,Any} = opts_regional_comparison(:blue, :cyan)
    opts_calibration::Dict{Symbol,Any} = opts_regional_comparison(:red, :magenta)
    f = Figure(size=fig_size)

    r_keys_validation = keys(regional_stats_validation)
    r_keys_calibration = keys(regional_stats_calibration)
    r_keys = union(r_keys_validation, r_keys_calibration)
    n_plots = length(r_keys)

    for (idx, key) in enumerate(r_keys)
        ax_title_validation = if key ∈ r_keys_validation
            "$(regional_stats_validation[key].n_locations) valid."
        else
            ""
        end

        ax_title_calibration = if key ∈ r_keys_calibration
            "$(regional_stats_calibration[key].n_locations) calib."
        else
            ""
        end

        ax_title_reefs = if ax_title_validation == ""
            ax_title_calibration
        elseif ax_title_calibration == ""
            ax_title_validation
        else
            ax_title_validation * " and " * ax_title_calibration
        end
        ax_title = "$(string(key))\n$ax_title_reefs" * " reefs"

        ax = Axis(
            f[fig_coord(idx, n_plots)...],
            title=ax_title,
            xlabel="Year",
            ylabel="Relative cover",
            limits=axis_limits,
            xticks=START_YEAR:2:END_YEAR,
            xticklabelrotation=π / 4
        )

        if key ∈ r_keys_validation
            plot_regional_comparison!(
                ax,
                regional_stats_validation[key].model_stats,
                replace(regional_stats_validation[key].historic_stats, missing => NaN),
                opts_validation,
            )

            n_locations = regional_stats_validation[key].n_locations
            ax_title * "$n_locations validation Reefs"
        end

        if key ∈ r_keys_calibration
            plot_regional_comparison!(
                ax,
                regional_stats_calibration[key].model_stats,
                replace(regional_stats_calibration[key].historic_stats, missing => NaN),
                opts_calibration,
            )

            n_locations = regional_stats_calibration[key].n_locations
            ax_title * "; $n_locations calibration Reefs"
        end

        ax.title[] = ax_title
    end

    legend_regional_comparison!(f[end+1, :], opts_validation, opts_calibration)

    Label(f[0, :], fig_title; fontsize=20)

    return f
end

function opts_regional_comparison(
    line_base_color::Symbol, historic_base_color::Symbol; alpha=0.1, markersize=7
)::Dict{Symbol,Any}
    return Dict(
        :model_line_color => line_base_color,
        :model_band_color => (line_base_color, alpha),
        :historic_line_color => historic_base_color,
        :historic_band_color => (historic_base_color, alpha),
        :markersize => markersize
    )
end

function legend_regional_comparison!(
    fig::Union{GridPosition,GridLayout,Figure};
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::Nothing
    model_leg_element = [
        LineElement(color=get(opts, :model_line_color, COLORS[:model_line_color])),
        PolyElement(color=get(opts, :model_band_color, COLORS[:model_band_color]))
    ]
    observed_leg_element = [
        LineElement(color=get(opts, :historic_line_color, COLORS[:historic_line_color])),
        PolyElement(color=get(opts, :historic_band_color, COLORS[:historic_band_color]))
    ]
    leg_labels = ["Observed", "Model"]

    Legend(
        fig,
        [observed_leg_element, model_leg_element],
        leg_labels,
        orientation=:horizontal
    )

    return nothing
end
function legend_regional_comparison!(
    fig::Union{GridPosition,GridLayout,Figure},
    opts_validation::Dict{Symbol,Any},
    opts_calibration::Dict{Symbol,Any})::Nothing
    model_leg_element_valid = [
        LineElement(color=opts_validation[:model_line_color]),
        PolyElement(color=opts_validation[:model_band_color])
    ]
    historic_leg_element_valid = [
        LineElement(color=opts_validation[:historic_line_color]),
        PolyElement(color=opts_validation[:historic_band_color])
    ]
    model_leg_element_calib = [
        LineElement(color=opts_calibration[:model_line_color]),
        PolyElement(color=opts_calibration[:model_band_color])
    ]
    historic_leg_element_calib = [
        LineElement(color=opts_calibration[:historic_line_color]),
        PolyElement(color=opts_calibration[:historic_band_color])
    ]

    leg_labels = [
        "CoralBlox validation Reefs",
        "LTMP validation Reefs",
        "CoralBlox calibration Reefs",
        "LTMP calibration Reefs"]

    leg_elements = [
        model_leg_element_valid,
        historic_leg_element_valid,
        model_leg_element_calib,
        historic_leg_element_calib
    ]

    Legend(
        fig[end+1, :],
        leg_elements,
        leg_labels,
        orientation=:horizontal
    )

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
    ref_years=START_YEAR:END_YEAR,
    fig_title="",
    fig_size=(800, 800)
)::Figure
    s_rac = (dropdims(sum(model_results.raw, dims=2), dims=2) .* site_k_area(dom)') ./ loc_area(dom)'

    north_res = ADRIA.DataCube(s_rac[:, region_masks[1]];
        timesteps=ref_years, locations=1:count(region_masks[1]))
    central_res = ADRIA.DataCube(s_rac[:, region_masks[2]];
        timesteps=ref_years, locations=1:count(region_masks[2]))
    south_res = ADRIA.DataCube(s_rac[:, region_masks[3]];
        timesteps=ref_years, locations=1:count(region_masks[3]))

    f = plot_all_regions(
        north_res, ltmp_north, central_res, ltmp_central, south_res, ltmp_south;
        fig_title=fig_title, fig_size=fig_size
    )

    return f
end
function plot_all_regions(
    north_results, north_obs, central_results, central_obs, south_results, south_obs;
    fig_title="", fig_size=(800, 500)
)::Figure
    f = Figure(; size=fig_size)

    plot_region!(
        f[1, 1],
        "North GBR\n$(size(north_results, 2)) validation locations",
        north_obs,
        north_results
    )
    plot_region!(
        f[1, 2],
        "Central GBR\n$(size(central_results, 2)) validation locations",
        central_obs,
        central_results
    )
    plot_region!(
        f[1, 3],
        "South GBR\n$(size(south_results, 2)) validation locations",
        south_obs,
        south_results
    )

    obs_el = [
        PolyElement(color=(:black, 0.4)),
        LineElement(color=:black)
    ]
    sim_el = [
        PolyElement(color=(:red, 0.4)),
        LineElement(color=:red)
    ]

    Legend(
        f[2, :],
        [obs_el, sim_el],
        ["LTMP", "CoralBlox"],
        orientation=:horizontal
    )

    # linkyaxes!(north_ax, central_ax, south_ax)
    # Label(f[0, :], fig_title; fontsize=20)

    # resize_to_layout!(f)
    return f
end

"""
    plot_region(f::Figure, row::Int64, col::Int64, title::String, obs::DataFrame, sim::YAXArray; showlegend::Bool = false, legend_row::Int64 = 2, legend_col::Int64 = 2)::Axis

"""
function plot_region!(
    g::GridPosition,
    title::String,
    obs::DataFrame,
    sim::YAXArray;
)::Nothing
    ax = Axis(
        g,
        title=title,
        xlabel="year",
        ylabel="Relative Absolute Cover",
        limits=(nothing, (0.0, 1.0))
    )

    # mean_agg = dropdims(mean(sim, dims=:locations), dims=:locations)
    confints = ADRIA.analysis.series_confint(read(sim))
    xs::Vector{Float64} = collect(sim.timesteps)

    band!(ax, obs.Year, obs.lower, obs.upper, color=(:black, 0.4))
    band!(ax, xs, confints[:, 1], confints[:, 3], color=(:red, 0.4))

    lines!(ax, obs.Year, obs.response, color=:black)
    lines!(ax, xs, confints[:, 2], color=:red, linewidth=2)

    return nothing
end
