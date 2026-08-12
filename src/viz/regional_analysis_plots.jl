"""
    save_regional_analysis_plots(raw_data, dom, calibration_store, test_store, out_dir,
        ltmp_north, ltmp_central, ltmp_south, north_mask, central_mask, south_mask)

Generate and save all regional comparison plots under `out_dir/regional_analysis/`.
"""
function save_regional_analysis_plots(
    raw_data::Array{Float64,4},
    dom,
    calibration_store::LocationDataStore,
    test_store::LocationDataStore,
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

    region_stats_test = region_stats(
        regions, region_masks, raw_data, dom; observations=test_store
    )
    region_stats_calibration = region_stats(
        regions, region_masks, raw_data, dom; observations=calibration_store
    )

    fig_regions_test = plot_regional_comparison(
        region_stats_test,
        keys(region_stats_test);
        opts=merge(base_opts, Dict{Symbol,Any}(:title => "Regional comparison - Test")),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 400))
    )
    save(joinpath(regional_analysis_save_dir, "regional__regions__test.png"), fig_regions_test)

    fig_regions_calib = plot_regional_comparison(
        region_stats_calibration,
        keys(region_stats_calibration);
        opts=merge(base_opts, Dict{Symbol,Any}(:title => "Regional comparison - Calibration")),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 400))
    )
    save(joinpath(regional_analysis_save_dir, "regional__regions__calibration.png"), fig_regions_calib)

    mareas = unique(canonical_reefs.management_area_short)
    marea_masks = [canonical_reefs.management_area_short .== area for area in mareas]
    r_stats_test = region_stats(
        Symbol.(mareas), marea_masks, raw_data, dom; observations=test_store
    )
    r_stats_calibration = region_stats(
        Symbol.(mareas), marea_masks, raw_data, dom; observations=calibration_store
    )

    fig_mang_area = plot_regional_comparison(
        r_stats_test,
        keys(r_stats_test);
        opts=merge(
            base_opts,
            Dict{Symbol,Any}(:title => "Management areas comparison\nTest Locations")
        ),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 400))
    )
    save(joinpath(regional_analysis_save_dir, "regional__management_areas__test.png"), fig_mang_area)

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
    spatial_group_stats_test = region_stats(
        Symbol.("Group " .* string.(spatial_groups)), spatial_group_masks, raw_data, dom;
        observations=test_store
    )

    keys_sort = sortperm(collect(parse.(Int, getindex.(split.(string.(keys(spatial_group_stats_test)), " "), 2))))

    fig_spat_gps = plot_regional_comparison(
        spatial_group_stats_test,
        keys(spatial_group_stats_test)[keys_sort];
        opts=merge(base_opts, Dict{Symbol,Any}(:title => "Reef groups (test reefs)", :invert_positions => true)),
        fig_opts=Dict{Symbol,Any}(:size => (1000, 800))
    )
    save(joinpath(regional_analysis_save_dir, "regional__spatial_grouping__test.png"), fig_spat_gps)

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

    loc_ids = test_store.domain_gpkg.UNIQUE_ID
    test_mask = (loc_ids .∈ Ref(test_store.ltmp_unique_ids))
    north_test_mask = north_mask .&& test_mask
    central_test_mask = central_mask .&& test_mask
    south_test_mask = south_mask .&& test_mask

    f_all_regions_test = plot_all_regions(
        dom, raw_data, ltmp_north, ltmp_central, ltmp_south;
        region_masks=[north_test_mask, central_test_mask, south_test_mask],
        fig_title="Regional comparison (modelled LTMP data)\nTest reefs",
        fig_size=(900, 400)
    )
    save(joinpath(regional_analysis_save_dir, "locs_reg_test.png"), f_all_regions_test)

    return nothing
end
