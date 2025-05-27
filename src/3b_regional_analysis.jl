include("analysis/regional_analysis.jl")

f_all_regions = plot_all_regions(
    dom, rs_raw
)
save(joinpath(OUT_DIR, "locs_reg.png"), f_all_regions)

# Plot all regions with validation locations only
loc_ids = VALIDATION_STORE.domain_gpkg.UNIQUE_ID
validation_mask = (loc_ids .∈ Ref(VALIDATION_STORE.ltmp_unique_ids))
NORTH_VALIDATION_MASK = NORTH_MASK .&& validation_mask
CENTRAL_VALIDATION_MASK = CENTRAL_MASK .&& validation_mask
SOUTH_VALIDATION_MASK = SOUTH_MASK .&& validation_mask

f_all_regions_validation = plot_all_regions(
    dom, rs_raw;
    region_masks=[NORTH_VALIDATION_MASK, CENTRAL_VALIDATION_MASK, SOUTH_VALIDATION_MASK]
)
save(joinpath(OUT_DIR, "locs_reg_validation.png"), f_all_regions_validation)

canonical_reefs = dom.loc_data

# * Managemant Areas
mareas = unique(canonical_reefs.management_area_short)
marea_masks = [canonical_reefs.management_area_short .== area for area in mareas]
r_stats = region_stats(Symbol.(mareas), marea_masks; observations=VALIDATION_STORE)

fig_mang_area = plot_regional_comparison(r_stats; fig_title="Management areas comparison\nValidation Locations", fig_size=(1000, 400))
save(joinpath(OUT_DIR, "regional__management_areas.png"), fig_mang_area)

# * Spatial Grouping
spatial_groups = unique(canonical_reefs.SPATIAL_GROUPING)
spatial_group_masks = [canonical_reefs.SPATIAL_GROUPING .== group for group in spatial_groups]
spatial_group_stats = region_stats(Symbol.("Group " .* string.(spatial_groups)), spatial_group_masks; observations=VALIDATION_STORE)

fig_spat_gps = plot_regional_comparison(spatial_group_stats; fig_title="Spatial grouping comparison\nValidation Locations", fig_size=(1000, 800))
save(joinpath(OUT_DIR, "regional__spatial_grouping.png"), fig_spat_gps)

# * Sectors
sectors_path = "C:/Users/pribeiro/AIMS/Datasets/AIMS_Sectors.csv"
sectors_df = CSV.read(sectors_path, DataFrame)

sector_names = collect(skipmissing(unique(sectors_df.aimssector)))
n_sectors = length(sector_names)
n_locations = size(canonical_reefs, 1)
sector_masks = fill(falses(n_locations), n_sectors)

for (i, name) in enumerate(sector_names)
    sector_ids = sectors_df.ReefID[BitVector(replace(sectors_df.aimssector .== name, missing => false))]
    sector_masks[i] = canonical_reefs.GBRMPA_ID .∈ Ref(sector_ids)
end

# * Bioregions

sectors_group_stats = region_stats(Symbol.(sector_names), sector_masks; observations=VALIDATION_STORE)
fig_sectors = plot_regional_comparison(sectors_group_stats; fig_title="Sectors comparison\nValidation Locations", fig_size=(1000, 600))
save(joinpath(OUT_DIR, "regional__sectors.png"), fig_sectors)
