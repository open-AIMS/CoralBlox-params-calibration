# * Reef groups scale
# * Spatial Grouping
include("../analysis/regional_analysis.jl")
spatial_groups = unique(domain_gpkg.CB_CALIB_GROUPS)
spatial_group_masks = [domain_gpkg.CB_CALIB_GROUPS .== g for g in spatial_groups]

spatial_group_stats_valid = region_stats(Symbol.("Group " .* string.(spatial_groups)), spatial_group_masks; observations=VALIDATION_STORE)
spatial_group_stats_calib = region_stats(Symbol.("Group " .* string.(spatial_groups)), spatial_group_masks; observations=CALIBRATION_STORE)

# Plot just valid. or calib.
data = [spatial_group_stats_valid, spatial_group_stats_calib]
reef_group_type = ["validation", "calibration"]

xticks_labels = string.(START_YEAR:END_YEAR)
year_labels = string.(Base.union(collect(START_YEAR:5:END_YEAR), [END_YEAR]))
xticks_labels[xticks_labels.∉Ref(year_labels)] .= ""
xticks = (START_YEAR:END_YEAR, xticks_labels)

yticknumbers = collect(0.0:0.2:1.0)
yticktext = getindex.(split.(string.(round.(yticknumbers .* 100)), "."), 1) .* "%"
ytick_labels = (yticknumbers, yticktext)

GeoMakie.GO.simplify.(geometries, ratio=0.5)

include("../plot/regional_comparison.jl")
for (idx_d, d) in enumerate(data)
    # idx_d = 1
    # d = data[idx_d]
    fig_reef_groups = plot_regional_comparison(
        d;
        opts=Dict{Symbol,Any}(
            :show_title => false,
            :titlesize => 11pt,
            :textlabelsize => 9pt,
            :textlabelbackground => "#ededeb",
            :showtextlabel => true
        ),
        fig_opts=Dict{Symbol,Any}(:size => (750, 750)),
        axis_opts=Dict{Symbol,Any}(
            :titlesize => 9pt,
            :xlabelsize => 9pt,
            :ylabelsize => 9pt,
            :xticklabelsize => 9pt,
            :yticklabelsize => 9pt,
            :xticklabelrotation => 0,
            :xticks => xticks,
            :yticks => ytick_labels,
        ),
        legend_opts=Dict{Symbol,Any}(:labelsize => 9pt)
    )

    save(fig_path * "/reef_groups_$(reef_group_type[idx_d]).png", fig_reef_groups; px_per_unit=(300 / inch))
end