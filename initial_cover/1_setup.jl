include("common.jl")

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

# Avoid reloading the domain every time
# Load ReefModDomain
if (!@isdefined(dom) || reload_domain)
    if reefmod_domain
        if start_year < 2008
            start_year = 2008
            @warn "Setting start year to $(start_year). 2008 is the earliest possible start for ReefModDomain."
        end

        @info "Loading ReefModDomain"
        dom = ADRIA.load_domain(ReefModDomain, reefmod_domain_path, "45", timeframe=(start_year, end_year))
    elseif !reefmod_domain
        # Load the RME Domain
        if start_year < 2000
            start_year = 2000
            @warn "Setting start year to $(start_year). 2000 is earlier possible start for RMEDomain."
        end
        @info "Loading RMEDomain"
        dom = ADRIA.load_domain(RMEDomain, rme_domain_path, "45", timeframe=(start_year, end_year))
    end

    @info "Attaching historic DHW"
    historic_dhw_path = joinpath(rme_domain_path, "data_files", "dhw", "GBR_past_DHW_CRW_5km_1985_2022_Dec_2022.csv")
    dhw_data_df = CSV.read(historic_dhw_path, DataFrame)

    # Available DHW data starts 1985 - 2022
    target_years = string.(start_year:end_year)
    dhw_data = reshape(Matrix(dhw_data_df[:, target_years])', 15, 3806, 1)

    dom.dhw_scens = ADRIA.DataCube(dhw_data; timesteps=target_years, locs=collect(caxes(dom.dhw_scens)[2]), scenarios=1:1)

    @info "Loading default parameters"
    scens = ADRIA.param_table(dom)

    reload_domain = false

    @info "Forcing Results Rerun"
    rerun = true
end

if !@isdefined(LTMP_DATA) || reload_ltmp
    LTMP_DATA = CSV.read("ltmp_data/modelled_brms.beta.ry.disp.csv", DataFrame, header=true)
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

    reload_tmp = false
end

if !@isdefined(region_shps) || reload_shp
    region_shps = GDF.read(ltmp_shp)
    reload_shp = false
end

if !@isdefined(NORTH_MASK)
    const NORTH_MASK = BitVector([AG.contains(region_shps.geometry[1], AG.centroid(polygn)) for polygn in dom.site_data.geom])  # .&& ltmp_loc_mask
    const CENTRAL_MASK = BitVector([AG.contains(region_shps.geometry[2], AG.centroid(polygn)) for polygn in dom.site_data.geom])  # .&& ltmp_loc_mask
    const SOUTH_MASK = BitVector([AG.contains(region_shps.geometry[3], AG.centroid(polygn)) for polygn in dom.site_data.geom])  # .&& ltmp_loc_mask
    const NOT_CONTAINED = (!).(NORTH_MASK .|| CENTRAL_MASK .|| SOUTH_MASK)
end


if !isdefined(Main, :uniform_initial_cover)
    uniform_initial_cover = false
end

if uniform_initial_cover
    @info "Using LTMP data to intialise uniform cover for north, central and southern GBR regions"
    north_index = findfirst(x -> x >= start_year, ltmp_north.Year)
    central_index = findfirst(x -> x >= start_year, ltmp_central.Year)
    south_index = findfirst(x -> x >= start_year, ltmp_south.Year)

    north_cover::Float64 = ltmp_north.response[north_index]
    central_cover::Float64 = ltmp_central.response[central_index]
    south_cover::Float64 = ltmp_south.response[south_index]

    interp_cover = (north_cover + central_cover + south_cover) / 3

    # Maintain species distributions, normalise to meet equivalent cover
    init_loc_cover = sum(dom.init_coral_cover, dims=:species)

    dom.init_coral_cover[locs=NORTH_MASK] .*= north_cover ./ init_loc_cover[locs=NORTH_MASK]
    dom.init_coral_cover[locs=CENTRAL_MASK] .*= central_cover ./ init_loc_cover[locs=CENTRAL_MASK]
    dom.init_coral_cover[locs=SOUTH_MASK] .*= south_cover ./ init_loc_cover[locs=SOUTH_MASK]

    # Locations not contained in the shapes defined, are assigned averaged values
    if any(NOT_CONTAINED)
        dom.init_coral_cover[locs=NOT_CONTAINED] .*= interp_cover ./ init_loc_cover[locs=NOT_CONTAINED]
    end

    dom.init_coral_cover ./= dom.site_data.k'
    # Convert total cover to relative cover

    uniform_initial_cover = false
end

if !@isdefined(north_res) && @isdefined(s_rac)
    north_res = s_rac[sites=NORTH_MASK]
    central_res = s_rac[sites=CENTRAL_MASK]
    south_res = s_rac[sites=SOUTH_MASK]
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

ltmp_reefmod_idxs = [
    ismissing(id) ? -1 : findfirst(
        x -> x == id,
        dom.site_data.UNIQUE_ID
    ) for id in ltmp_reef_data.RME_UNIQUE_ID
]

ENV["ADRIA_DEBUG"] = false
