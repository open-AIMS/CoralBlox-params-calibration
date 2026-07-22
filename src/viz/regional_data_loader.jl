"""
    RegionalAnalysisData

Analysis-only LTMP regional observations and location region masks, consumed by
`perf_metrics.jl`'s regional stats and the regional/location-timeseries plots in this
module. Loaded exclusively via [`load_regional_analysis_data`](@ref), never by the
calibration entry point (`run_calibration` never touches regional splits).
"""
struct RegionalAnalysisData
    ltmp_north::DataFrame
    ltmp_central::DataFrame
    ltmp_south::DataFrame
    north_mask::BitVector
    central_mask::BitVector
    south_mask::BitVector
    not_contained::BitVector
end

function _ltmp_region_data(region_name::String, ltmp_data::DataFrame)::DataFrame
    region_mask = [region_name == reg for reg in ltmp_data.Region]
    return ltmp_data[region_mask, :]
end

function _region_shape_mask(dom, region_shapes, idx::Int)::BitVector
    region_shape = region_shapes.geometry[idx]
    geoms = GI.geometry.(eachrow(dom.loc_data))
    return [AG.contains(region_shape, AG.centroid(geom)) for geom in geoms]
end

"""
    load_regional_analysis_data(dom, ltmp_modelled_obs_path, ltmp_shp_path)::RegionalAnalysisData

Load LTMP modelled observations split by region ("Northern GBR"/"Central GBR"/"Southern
GBR") and compute per-location region containment masks against `dom.loc_data`.
"""
function load_regional_analysis_data(
    dom, ltmp_modelled_obs_path::String, ltmp_shp_path::String
)::RegionalAnalysisData
    ltmp_data = CSV.read(ltmp_modelled_obs_path, DataFrame, header=true)
    ltmp_data[!, :Region] = String.(ltmp_data[:, :Region])
    regions = ["Northern GBR", "Central GBR", "Southern GBR"]
    ltmp_north, ltmp_central, ltmp_south = _ltmp_region_data.(regions, [ltmp_data])

    region_shps = GDF.read(ltmp_shp_path)
    north_mask = _region_shape_mask(dom, region_shps, 1)
    central_mask = _region_shape_mask(dom, region_shps, 2)
    south_mask = _region_shape_mask(dom, region_shps, 3)
    not_contained = (!).(north_mask .|| central_mask .|| south_mask)

    return RegionalAnalysisData(
        ltmp_north, ltmp_central, ltmp_south,
        north_mask, central_mask, south_mask, not_contained
    )
end
