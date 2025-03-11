include("./common/common.jl")
include("./common/cover_construction.jl")

if !@isdefined(canonical_gpkg)
    @info "Loading Canonical gpkg"
    canonical_gpkg = GDF.read(CANONICAL_PATH)
    ltmp_loc_mask = canonical_gpkg.is_LTMP_reef .!= 0
end

# Avoid reloading the domain every time
# Load ReefModDomain
if (!@isdefined(dom) || reload_domain)
    if START_YEAR < 2008
        START_YEAR = 2008
        @warn "Setting start year to $(START_YEAR). 2008 is the earliest possible start for ReefModDomain."
    end

    @info "Loading ReefModDomain"
    dom = ADRIA.load_domain(ReefModDomain, REEFMOD_DOMAIN_PATH, "45", timeframe=(START_YEAR, END_YEAR))

    @info "Attaching historic DHW"
    dhw_data_df = CSV.read(HISTORIC_DHW_PATH, DataFrame)

    # Available DHW data starts 1985 - 2022
    target_years = string.(START_YEAR:END_YEAR)
    locs = collect(dom.dhw_scens.locs)

    n_timesteps = length(target_years)
    n_locs = length(locs)

    dhw_data = reshape(Matrix(dhw_data_df[:, target_years])', n_timesteps, n_locs, 1)
    dom.dhw_scens = ADRIA.DataCube(dhw_data; timesteps=target_years, locs=locs, scenarios=1:1)

    new_cyclone_mortality_scens = open_dataset("../ltmp_data/new_cyclone_mortality_scens.nc").cyclone_mortality_scens
    dom.cyclone_mortality_scens .= new_cyclone_mortality_scens

    reload_domain = false
end

function ltmp_period(ltmp_region_name::String, ltmp_data::DataFrame, START_YEAR::Int64, END_YEAR::Int64)::BitVector
    region_tf = ltmp_data[ltmp_data.Region.==ltmp_region_name, :Year]
    return (region_tf .>= START_YEAR) .& (region_tf .<= END_YEAR)
end

function ltmp_modelled_results(ltmp_region_name::String, ltmp_data::DataFrame)::DataFrame
    region_mask = [ltmp_region_name == reg for reg in ltmp_data.Region]
    return ltmp_data[region_mask, :]
end

if !@isdefined(ltmp_data)
    ltmp_data = CSV.read(LTMP_MODELLED_OBS_PATH, DataFrame, header=true)
    ltmp_data[!, :Region] = String.(ltmp_data[:, :Region])
    regions = ["Northern GBR", "Central GBR", "Southern GBR"]

    ltmp_north, ltmp_central, ltmp_south = ltmp_modelled_results.(
        regions, [ltmp_data]
    )

    ltmp_north_period, ltmp_central_period, ltmp_south_period =
        ltmp_period.(regions, [ltmp_data], [START_YEAR], [END_YEAR])
end

if !@isdefined(region_shps)
    region_shps = GDF.read(LTMP_SHP_PATH)
end

function _region_shape_mask(dom, region_shapes, idx)::BitVector
    region_shape = region_shapes.geometry[idx]
    return [AG.contains(region_shape, AG.centroid(poly)) for poly in dom.loc_data.geom]
end

if !@isdefined(NORTH_MASK)
    const NORTH_MASK = _region_shape_mask(dom, region_shps, 1)  # .&& ltmp_loc_mask
    const CENTRAL_MASK = _region_shape_mask(dom, region_shps, 2)  # .&& ltmp_loc_mask
    const SOUTH_MASK = _region_shape_mask(dom, region_shps, 3)  # .&& ltmp_loc_mask
    const NOT_CONTAINED = (!).(NORTH_MASK .|| CENTRAL_MASK .|| SOUTH_MASK)
end

location_classification = CSV.read(LOC_CLASS_PATH, DataFrame)

# ? Delete
# ? n_classifications = maximum(location_classification.consecutive_classification)

# set ADRIA Env variable to prevent an error during run_model
ENV["ADRIA_DEBUG"] = false
