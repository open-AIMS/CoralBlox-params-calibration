"""
# Example
```
management_areas = unique(dom.loc_data.management_area_short)
management_area_masks = [dom.loc_data.management_area_short .== area for area in management_areas]

management_area_stats = region_stats(Symbol.(management_area_masks), management_areas)

plot_regional_comparison(management_area_stats)
```
"""
function plot_regional_comparison(
    regional_stats::NamedTuple;
    fig_size=(800, 800), axis_limits=(nothing, (0.0, 1.0)), fig_title=""
)
    r_keys = keys(regional_stats)

    band_alpha = 0.2
    markersize = 7

    n_plots = length(r_keys)
    # ntimesteps = size(regional_stats[r_keys[1]].model_stats, 2)
    xdata = START_YEAR:END_YEAR

    f = Figure(size=fig_size)
    for (idx, key) in enumerate(r_keys)
        model_stats = regional_stats[key].model_stats
        historic_stats = replace(regional_stats[key].historic_stats, missing => NaN)
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

        scatterlines!(ax, xdata, model_stats[2, :], color=:blue, markersize=markersize)
        band!(ax, xdata, model_stats[1, :], model_stats[3, :], color=(:blue, band_alpha))

        scatterlines!(ax, xdata, historic_stats[2, :], color=:red, markersize=markersize)
        band!(ax, xdata, historic_stats[1, :], historic_stats[3, :]; color=(:red, band_alpha))
    end

    model_leg_element = [LineElement(color=:blue), PolyElement(color=(:blue, band_alpha))]
    historic_leg_element = [LineElement(color=:red), PolyElement(color=(:red, band_alpha))]
    leg_labels = ["Model data", "Historic Data"]
    Legend(f[end+1, :], [model_leg_element, historic_leg_element], leg_labels, orientation=:horizontal)

    Label(f[0, :], fig_title; fontsize=20)

    return f
end
