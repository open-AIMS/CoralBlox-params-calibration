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
    if start_year < 2008
        start_year = 2008
        @warn "Setting start year to $(start_year). 2008 is the earliest possible start for ReefModDomain."
    end

    @info "Loading ReefModDomain"
    dom = ADRIA.load_domain(ReefModDomain, REEFMOD_DOMAIN_PATH, "45", timeframe=(start_year, end_year))

    @info "Attaching historic DHW"
    dhw_data_df = CSV.read(HISTORIC_DHW_PATH, DataFrame)

    # Available DHW data starts 1985 - 2022
    target_years = string.(start_year:end_year)
    locs = collect(caxes(dom.dhw_scens)[2])

    n_timesteps = length(target_years)
    n_locs = length(locs)

    dhw_data = reshape(Matrix(dhw_data_df[:, target_years])', n_timesteps, n_locs, 1)
    dom.dhw_scens = ADRIA.DataCube(dhw_data; timesteps=target_years, locs=locs, scenarios=1:1)

    @info "Loading default parameters"
    scens = ADRIA.param_table(dom)

    reload_domain = false
end

function ltmp_period(ltmp_region_name::String, ltmp_data::DataFrame, start_year::Int64, end_year::Int64)::BitVector
    region_tf = ltmp_data[ltmp_data.Region .== ltmp_region_name, :Year]
    return (region_tf .>= start_year) .& (region_tf .<= end_year)
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
        ltmp_period.(regions, [ltmp_data], [start_year], [end_year])
end

if !@isdefined(region_shps)
    region_shps = GDF.read(LTMP_SHP_PATH)
end

if !@isdefined(NORTH_MASK)
    const NORTH_MASK = BitVector([AG.contains(region_shps.geometry[1], AG.centroid(polygn)) for polygn in dom.loc_data.geom])  # .&& ltmp_loc_mask
    const CENTRAL_MASK = BitVector([AG.contains(region_shps.geometry[2], AG.centroid(polygn)) for polygn in dom.loc_data.geom])  # .&& ltmp_loc_mask
    const SOUTH_MASK = BitVector([AG.contains(region_shps.geometry[3], AG.centroid(polygn)) for polygn in dom.loc_data.geom])  # .&& ltmp_loc_mask
    const NOT_CONTAINED = (!).(NORTH_MASK .|| CENTRAL_MASK .|| SOUTH_MASK)
end

if !@isdefined(north_res) && @isdefined(s_rac)
    north_res = s_rac[locs=NORTH_MASK]
    central_res = s_rac[locs=CENTRAL_MASK]
    south_res = s_rac[locs=SOUTH_MASK]
end

location_classification = CSV.read(LOC_CLASS_PATH, DataFrame)
n_classifications = maximum(location_classification.consecutive_classification)

# Load manta observations for reef location classes
manta_tow_classes = open_dataset(LOC_CLASS_TARGET_PATH)

# Force memory load
manta_tow_mean = readcubedata(manta_tow_classes.mean)
manta_tow_std = readcubedata(manta_tow_classes.std)

# Load manta tow ltmp reef level data
ltmp_reef_data = GDF.read(LTMP_REEF_DATA_PATH)

# Order year columns in ascending order
ltmp_reef_years = parse.(Int64, names(ltmp_reef_data)[5:end])
ltmp_reef_perm = sortperm(ltmp_reef_years) .+ 4

ltmp_reef_data_names = names(ltmp_reef_data)
ltmp_reef_data_names[5:end] .= ltmp_reef_data_names[ltmp_reef_perm]

# Rorder columns
ltmp_reef_data = select!(ltmp_reef_data, ltmp_reef_data_names...)

# Rescale to be proportions
ltmp_reef_data[:, 5:end] ./= 100
first_yr_idx = findfirst(x -> x == "2008", names(ltmp_reef_data))
raw_ltmp_reef_data = Matrix(ltmp_reef_data[:, first_yr_idx:end])

# For each ltmp location calculate the row index for the domain
all_ltmp_idxs = [
    ismissing(id) ? -1 : findfirst(
        x -> x == id,
        dom.loc_data.UNIQUE_ID
    ) for id in ltmp_reef_data.RME_UNIQUE_ID
]
all_ltmp_reef = copy(raw_ltmp_reef_data)

# Calibration Locations
limited_locations = ["16015100104", "16025100104", "14114100104", "18075100104"]
location_names = ["Mackay Reef", "Opal Reef", "Macgillivray Reef", "John Brewer Reef"]
limited_loc_idxs = [findfirst(x -> !ismissing(x) && x == id, ltmp_reef_data.RME_UNIQUE_ID) for id in limited_locations]
raw_ltmp_reef_data = raw_ltmp_reef_data[limited_loc_idxs, :]

composition_data = open_dataset(COMPOSITION_PATH)

# For each target location get its row index in the domain
target_dom_idxs = [findfirst(x -> x == id, dom.loc_data.UNIQUE_ID) for id in limited_locations]

# Extract the ltmp data for each target location
temporal_range = start_year:end_year
# [year ⋅ taxa ⋅ locs]
rm_ltmp_taxa = Array{Union{Missing,Float64}}(missing, length(temporal_range), 5, 4)
for (j, loc_name) in enumerate(location_names)
    rm_ltmp_taxa[:, :, j] .= composition_data.mean[location=At(loc_name)].data[(2008-1992+1):(2008-1992)+length(temporal_range), :]
end

dom.loc_data[!, :depth_med] .= canonical_gpkg.depth_med

# set ADRIA Env variable to prevent an error during run_model
ENV["ADRIA_DEBUG"] = false
