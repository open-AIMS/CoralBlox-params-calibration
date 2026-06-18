include("common.jl")

HISTORICAL_DHW_PATH = joinpath(
    @__DIR__, "../src/historical_dhw_data_gen/data/historical_dhw.nc"
)
HISTORICAL_CYCLONE_MORTALITY_PATH = joinpath(
    @__DIR__,
    "../src/historical_disturbance_mortality_data_gen/data/",
    "historical_disturbance_mortality_rates/historical_disturbance_mortality_rates.nc"
)

if !@isdefined(canonical_gpkg) || reload_canonical
    @info "Loading Canonical gpkg"
    canonical_gpkg = GDF.read(canonical_path)
    ltmp_loc_mask = canonical_gpkg.is_LTMP_reef .!= 0
    reload_canonical = false
end

if !@isdefined(reload_ltmp)
    reload_ltmp = true
end

if !@isdefined(reload_shp)
    reload_shp = true
end

if (!@isdefined(dom) || reload_domain)
    @info "Loading RMEDomain"
    dom = ADRIA.load_domain(RMEDomain, rme_domain_path, "45", timeframe=(start_year, end_year))

    @info "Attaching historic DHW and Cyclone/COTS data"
    new_dhw_scens = open_dataset(HISTORICAL_DHW_PATH).dhw_scens
    dom.dhw_scens .= read(new_dhw_scens[timesteps=At(start_year:end_year), scenarios=1])

    new_cyclone_mortality_scens = open_dataset(HISTORICAL_CYCLONE_MORTALITY_PATH).disturbance_mortality_scens
    dom.cyclone_mortality_scens .= read(new_cyclone_mortality_scens[:, :, :, [1]])

    reload_domain = false
    rerun = true
end

if !@isdefined(LTMP_DATA) || reload_ltmp
    LTMP_DATA = CSV.read(joinpath(@__DIR__, "../datasets/ltmp_data/modelled_brms.beta.ry.disp.csv"), DataFrame, header=true)
    LTMP_DATA[!, :Region] = String.(LTMP_DATA[:, :Region])

    ltmp_north_mask = ["Northern GBR" == reg for reg in LTMP_DATA.Region]
    ltmp_central_mask = ["Central GBR" == reg for reg in LTMP_DATA.Region]
    ltmp_south_mask = ["Southern GBR" == reg for reg in LTMP_DATA.Region]

    ltmp_north = LTMP_DATA[ltmp_north_mask, :]
    ltmp_central = LTMP_DATA[ltmp_central_mask, :]
    ltmp_south = LTMP_DATA[ltmp_south_mask, :]

    ltmp_north_period = (ltmp_north.Year .>= start_year) .& (ltmp_north.Year .<= end_year)
    ltmp_central_period = (ltmp_central.Year .>= start_year) .& (ltmp_central.Year .<= end_year)
    ltmp_south_period = (ltmp_south.Year .>= start_year) .& (ltmp_south.Year .<= end_year)

    reload_ltmp = false
end

if !@isdefined(region_shps) || reload_shp
    region_shps = GDF.read(ltmp_shp)
    reload_shp = false
end

function _region_shape_mask(dom, region_shapes, idx)::BitVector
    region_shape = region_shapes.geometry[idx]
    geoms = GI.geometry.(eachrow(dom.loc_data))
    return [AG.contains(region_shape, AG.centroid(geom)) for geom in geoms]
end

if !@isdefined(NORTH_MASK)
    const NORTH_MASK = _region_shape_mask(dom, region_shps, 1)
    const CENTRAL_MASK = _region_shape_mask(dom, region_shps, 2)
    const SOUTH_MASK = _region_shape_mask(dom, region_shps, 3)
    const NOT_CONTAINED = (!).(NORTH_MASK .|| CENTRAL_MASK .|| SOUTH_MASK)
end

location_classification = CSV.read(classification_path, DataFrame)
n_classifications = maximum(location_classification.consecutive_classification)

# Load manta observations for reef location classes
manta_tow_classes = open_dataset(manta_tow_class_path)

# Force memory load
manta_tow_mean = readcubedata(manta_tow_classes.mean)
manta_tow_std = readcubedata(manta_tow_classes.std)

# Load manta tow ltmp reef level data
ltmp_reef_data = GDF.read(ltmp_reef_data_path)

# Identify year columns by name (some cols like "match_reason" are not years)
let all_names = names(ltmp_reef_data)
    is_year = [tryparse(Int64, n) !== nothing for n in all_names]
    year_idxs = findall(is_year)
    non_year_idxs = findall(.!is_year)
    sorted_year_idxs = year_idxs[sortperm(parse.(Int64, all_names[year_idxs]))]
    global ltmp_reef_data = DataFrames.select(ltmp_reef_data, [non_year_idxs; sorted_year_idxs]...)
end
# Divide all year columns by 100 (stored as percentages)
first_yr_idx = findfirst(x -> tryparse(Int64, x) !== nothing, names(ltmp_reef_data))
ltmp_reef_data[:, first_yr_idx:end] ./= 100

first_yr_idx = findfirst(x -> x == "2008", names(ltmp_reef_data))
raw_ltmp_reef_data = Matrix(ltmp_reef_data[:, first_yr_idx:end])

ltmp_reefmod_idxs = [
    ismissing(id) ? -1 : findfirst(
        x -> x == id,
        dom.loc_data.UNIQUE_ID
    ) for id in ltmp_reef_data.RME_UNIQUE_ID
]

ENV["ADRIA_DEBUG"] = false

@info "Successfuly finished running 1_setup.jl"
