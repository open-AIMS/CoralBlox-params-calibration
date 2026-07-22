"""
    save_regional_analysis_plots(raw_data, dom, calibration_store, validation_store, out_dir,
        ltmp_north, ltmp_central, ltmp_south, north_mask, central_mask, south_mask)

Generate and save all regional comparison plots under `out_dir/regional_analysis/`.
"""
function save_regional_analysis_plots(
    raw_data::Array{Float64,4},
    dom,
    calibration_store::LocationDataStore,
    validation_store::LocationDataStore,
    out_dir::String,
    ltmp_north::DataFrame,
    ltmp_central::DataFrame,
    ltmp_south::DataFrame,
    north_mask::BitVector,
    central_mask::BitVector,
    south_mask::BitVector,
)::Nothing
    regional_analysis_save_dir = joinpath(out_dir, "regional_analysis")
    mkpath(regional_analysis_save_dir)

    canonical_reefs = dom.loc_data

    base_opts = Dict{Symbol,Any}(
        :textlabelbackground => "#ededeb",
        :showtextlabel => true,
        :invert_positions => false
    )

    region_masks = [north_mask, central_mask, south_mask]
    regions = [:North, :Central, :South]

    region_stats_validation = region_stats(
        regions, region_masks, raw_data, dom; observations=validation_store
    )
    region_stats_calibration = region_stats(
        regions, region_masks, raw_data, dom; observations=calibration_store
    )

    fig_regions_valid = plot_regional_comparison(
        region_stats_validation,
        keys(region_stats_validation);
        opts=merge(base_opts, Dict{Symbol,Any}(:title => "Regional comparison - Validation")),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 400))
    )
    save(joinpath(regional_analysis_save_dir, "regional__regions__validation.png"), fig_regions_valid)

    fig_regions_calib = plot_regional_comparison(
        region_stats_calibration,
        keys(region_stats_calibration);
        opts=merge(base_opts, Dict{Symbol,Any}(:title => "Regional comparison - Calibration")),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 400))
    )
    save(joinpath(regional_analysis_save_dir, "regional__regions__calibration.png"), fig_regions_calib)

    mareas = unique(canonical_reefs.management_area_short)
    marea_masks = [canonical_reefs.management_area_short .== area for area in mareas]
    r_stats_validation = region_stats(
        Symbol.(mareas), marea_masks, raw_data, dom; observations=validation_store
    )
    r_stats_calibration = region_stats(
        Symbol.(mareas), marea_masks, raw_data, dom; observations=calibration_store
    )

    fig_mang_area = plot_regional_comparison(
        r_stats_validation,
        keys(r_stats_validation);
        opts=merge(
            base_opts,
            Dict{Symbol,Any}(:title => "Management areas comparison\nValidation Locations")
        ),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 400))
    )
    save(joinpath(regional_analysis_save_dir, "regional__management_areas__validation.png"), fig_mang_area)

    fig_mang_area = plot_regional_comparison(
        r_stats_calibration,
        keys(r_stats_calibration);
        opts=merge(base_opts, Dict{Symbol,Any}(:title => "Management areas comparison\nCalibration Locations")),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 400))
    )
    save(joinpath(regional_analysis_save_dir, "regional__management_areas__calibration.png"), fig_mang_area)

    spatial_groups = unique(canonical_reefs.CB_CALIB_GROUPS)
    spatial_group_masks = [canonical_reefs.CB_CALIB_GROUPS .== group for group in spatial_groups]

    spatial_group_stats_calib = region_stats(
        Symbol.("Group " .* string.(spatial_groups)), spatial_group_masks, raw_data, dom;
        observations=calibration_store
    )
    spatial_group_stats_valid = region_stats(
        Symbol.("Group " .* string.(spatial_groups)), spatial_group_masks, raw_data, dom;
        observations=validation_store
    )

    keys_sort = sortperm(collect(parse.(Int, getindex.(split.(string.(keys(spatial_group_stats_valid)), " "), 2))))

    fig_spat_gps = plot_regional_comparison(
        spatial_group_stats_valid,
        keys(spatial_group_stats_valid)[keys_sort];
        opts=merge(base_opts, Dict{Symbol,Any}(:title => "Reef groups (validation reefs)", :invert_positions => true)),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 800))
    )
    save(joinpath(regional_analysis_save_dir, "regional__spatial_grouping__validation.png"), fig_spat_gps)

    fig_spat_gps = plot_regional_comparison(
        spatial_group_stats_calib,
        keys(spatial_group_stats_calib)[keys_sort];
        opts=merge(base_opts, Dict{Symbol,Any}(:title => "Reef groups (calibration reefs)", :invert_positions => true)),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 800))
    )
    save(joinpath(regional_analysis_save_dir, "regional__spatial_grouping__calibration.png"), fig_spat_gps)

    f_all_regions = plot_all_regions(
        dom, raw_data, ltmp_north, ltmp_central, ltmp_south;
        region_masks=[north_mask, central_mask, south_mask],
        fig_title="Regional comparison (modelled LTMP data)\nAll reefs",
        fig_size=(900, 400)
    )
    save(joinpath(regional_analysis_save_dir, "locs_reg.png"), f_all_regions)

    loc_ids = validation_store.domain_gpkg.UNIQUE_ID
    validation_mask = (loc_ids .∈ Ref(validation_store.ltmp_unique_ids))
    north_validation_mask = north_mask .&& validation_mask
    central_validation_mask = central_mask .&& validation_mask
    south_validation_mask = south_mask .&& validation_mask

    f_all_regions_validation = plot_all_regions(
        dom, raw_data, ltmp_north, ltmp_central, ltmp_south;
        region_masks=[north_validation_mask, central_validation_mask, south_validation_mask],
        fig_title="Regional comparison (modelled LTMP data)\nValidation reefs",
        fig_size=(900, 400)
    )
    save(joinpath(regional_analysis_save_dir, "locs_reg_validation.png"), f_all_regions_validation)

    return nothing
end
