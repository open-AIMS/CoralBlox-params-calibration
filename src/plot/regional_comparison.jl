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
    scatterlines!(
        ax, xdata, collect(model_stats[2, :]);
        color=opts[:model_line_color], markersize=opts[:markersize]
    )
    band!(ax, xdata, model_stats[1, :], model_stats[3, :]; color=opts[:model_band_color])

    scatterlines!(
        ax, xdata, historic_stats[2, :];
        color=opts[:historic_line_color], markersize=opts[:markersize]
    )
    band!(ax, xdata, historic_stats[1, :], historic_stats[3, :]; color=opts[:historic_band_color])

    return nothing
end
function plot_regional_comparison(
    regional_stats::NamedTuple;
    fig_size=(800, 800), axis_limits=(nothing, (0.0, 1.0)), fig_title=""
)
    opts::Dict{Symbol,Any} = opts_regional_comparison(:blue, :red)
    f = Figure(size=fig_size)
    r_keys = keys(regional_stats)
    n_plots = length(r_keys)

    for (idx, key) in enumerate(r_keys)
        n_locations = regional_stats[key].n_locations
        ax = Axis(
            f[fig_coord(idx, n_plots)...],
            title="$(string(key)) ($n_locations locations)",
            xlabel="Year",
            ylabel="Relative cover",
            limits=axis_limits,
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

    legend_regional_comparison!(f, opts)
    Label(f[0, :], fig_title; fontsize=20)

    return f
end

function plot_regional_comparison(
    regional_stats_validation::NamedTuple,
    regional_stats_calibration::NamedTuple;
    fig_size=(800, 800), axis_limits=(nothing, (0.0, 1.0)), fig_title=""
)
    opts_validation::Dict{Symbol,Any} = opts_regional_comparison(:blue, :cyan)
    opts_calibration::Dict{Symbol,Any} = opts_regional_comparison(:red, :magenta)
    f = Figure(size=fig_size)

    r_keys_validation = keys(regional_stats_validation)
    r_keys_calibration = keys(regional_stats_calibration)
    r_keys = union(r_keys_validation, r_keys_calibration)
    n_plots = length(r_keys)

    for (idx, key) in enumerate(r_keys)
        ax_title = "$(string(key))\n"
        ax = Axis(
            f[fig_coord(idx, n_plots)...],
            # title="$(string(key)) ($n_locations locations)",
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
            ax_title * "$n_locations validation locs"
        end

        if key ∈ r_keys_calibration
            plot_regional_comparison!(
                ax,
                regional_stats_calibration[key].model_stats,
                replace(regional_stats_calibration[key].historic_stats, missing => NaN),
                opts_calibration,
            )

            n_locations = regional_stats_calibration[key].n_locations
            ax_title * "; $n_locations calibration locs"
        end

        ax.title[] = ax_title
    end

    legend_regional_comparison!(f, opts_validation, opts_calibration)

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

function legend_regional_comparison!(fig::Figure, opts::Dict{Symbol,Any})::Nothing
    model_leg_element = [
        LineElement(color=opts[:model_line_color]),
        PolyElement(color=opts[:model_band_color])
    ]
    historic_leg_element = [
        LineElement(color=opts[:historic_line_color]),
        PolyElement(color=opts[:historic_band_color])
    ]
    leg_labels = ["Model data", "Historic Data"]

    Legend(
        fig[end+1, :],
        [model_leg_element, historic_leg_element],
        leg_labels,
        orientation=:horizontal
    )

    return nothing
end
function legend_regional_comparison!(
    fig::Figure,
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
        "Validation model data ",
        "Validation historic data",
        "Calibration model data",
        "Calibration historic data"]

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
