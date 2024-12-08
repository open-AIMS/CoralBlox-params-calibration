using StatsBase
using YAXArrays, DataFrames, NetCDF
using ADRIA: GDF

include("../common/perf_metrics.jl")

valid_historic_data_type() = (:ltmp, :transect)

"""
    rmse_scores(ltmp_data_path, rme_benchmark_path, canonical_path)

RMSE (Root Mean Square Error) for both ReefModEngine and LTMP average historic data
benchmark. In both cases, the error was computed against LTMP historic data. LTMP data for
locations for which the field RME_UNIQUE_ID was missing and also for locations with
duplicated data was discarded. Data for locations for which there was less than 4 years of
data collected was also discarded. For each location this selects the first year for which
there is LTMP data within the analyzed timeframe (2008-2022) and looks for the RME scenario
with the closest value for that year; selects only RME years for which there is LTMP data;
computes the RMSE score for RME data points and LTMP data; and computes the RMSE score for
LTMP benchmark (average of LTMP data for those years and that location) and LTMP data.

# Arguments
- `ltmp_data_path` : Path to manta tow LTMP geopackage (manta_tow_data_reef_lvl.gpkg)
- `rme_benchmark_path` : Path to RME model cover results from 2008 to 2022
- `canonical_path` : Path to canonical GBR geopackage (canonical_gbr_2024-04-23.gpkg)
- `ltmp_cover` : DataFrame with LTMP
- `rme_cover` : YAXArray with ReefModEngine model cover

# Returns
NamedTuple with three elements: rme_scores (vector of RMSE scores for RME model results and
LTMP historic data for each location), benchmark_scores (vector of RMSE scores for
benchmark, LTMP average, and LTMP historic data), n_observations (vector with number of
observations for each location)
"""
function rmse_scores(
    historic_data_path::String,
    reefmod_path::String,
    canonical_path::String,
    historic_data_type::Symbol
)::NamedTuple
    if historic_data_type ∉ valid_historic_data_type()
        ArgumentError("`historic_data_type` is $historic_data_type. "*
        "Valid options are $(valid_historic_data_type())")
    end

    canonical_gpkg = GDF.read(canonical_path)
    historic_data::YAXArray = if historic_data_type==:ltmp
        _ltmp_cube(historic_data_path)
    else
        _transect_cube(historic_data_path, canonical_gpkg)
    end
    reefmod_cover::YAXArray = _reefmod_cover(reefmod_path, historic_data, canonical_gpkg)

    return rmse_scores(historic_data, reefmod_cover)
end
function rmse_scores(historic_data::YAXArray, reefmod_cover::YAXArray)::NamedTuple
    target_timesteps::Vector{Int64} = _intersect_timerange(historic_data, reefmod_cover)
    historic_data = _filter_low_data_locations(historic_data[timesteps=At(target_timesteps)])

    n_locations = size(historic_data, :locations)
    reefmod_scores::Vector{Float64} = zeros(Float64, n_locations)
    benchmark_scores::Vector{Float64} = zeros(Float64, n_locations)
    n_observations::Vector{Int64} = zeros(Float64, n_locations)
    for (loc_id, ltmp_loc_cover) in enumerate(eachcol(historic_data))
        # Identify LTMP years with data
        years_with_data::BitVector = .!ismissing.(collect(ltmp_loc_cover))
        ltmp_target_cover = ltmp_loc_cover[years_with_data]
        target_years = target_timesteps[years_with_data]

        # Fill scores
        reefmod_scores[loc_id] = _reefmod_scores(reefmod_cover, ltmp_target_cover, target_years, loc_id)
        benchmark_scores[loc_id] = _benchmark_scores(ltmp_target_cover)
        n_observations[loc_id] = length(ltmp_target_cover)
    end
    return (
        reefmod_scores=reefmod_scores,
        benchmark_scores=benchmark_scores,
        n_observations=n_observations
    )
end

function _benchmark_scores(target_cover::Union{YAXArray, DataFrameRow})
    benchmark_cover = fill(mean(target_cover), length(target_cover))
    _benchmark_scores(benchmark_cover, target_cover)
end
function _benchmark_scores(benchmark_cover, target_cover::YAXArray)
    rmse(benchmark_cover, target_cover.data)
end

function _reefmod_scores(
    rme_cover::YAXArray,
    target_cover::Union{YAXArray, DataFrameRow},
    target_years::Vector{Int64},
    loc_idx::Int64
)::Float64
    rme_target_cover = rme_cover[timestep=At(target_years), location=loc_idx]
    rme_target_cover = _rme_closest_scenario(rme_target_cover, target_cover)
    rmse(rme_target_cover, target_cover.data)
end

function _intersect_timerange(historic_data::YAXArray, rme_cover::YAXArray)::Vector{Int64}
    start_year = max(minimum(historic_data.timesteps), minimum(rme_cover.timestep))
    end_year = min(maximum(historic_data.timesteps), maximum(rme_cover.timestep))
    return start_year:end_year
end

function _filter_low_data_locations(datacube::YAXArray; min_n_data::Int64=4)
    has_data = .!ismissing.(datacube)
    target_locs = dropdims(sum(has_data, dims=:timesteps) .>= min_n_data, dims=:timesteps)
    return datacube[locations=target_locs]
end

function _rme_closest_scenario(rme_target_cover::YAXArray, ltmp_loc_cover::DataFrameRow)::YAXArray
    first_year_with_ltmp_data = findmin(names(ltmp_loc_cover))[1]
    ltmp_first_year_cover = ltmp_loc_cover[first_year_with_ltmp_data]
    closest_scenario_idx = argmin(abs.(ltmp_first_year_cover .- rme_target_cover[1, :]))
    return rme_target_cover[:, closest_scenario_idx]
end
function _rme_closest_scenario(rme_target_cover::YAXArray, transect_target_cover::YAXArray)::YAXArray
    first_transect_year = minimum(transect_target_cover.timesteps)
    transect_first_year_cover = transect_target_cover[timesteps=At(first_transect_year)]
    closest_scenario_idx = argmin(abs.(transect_first_year_cover .- rme_target_cover[1, :]))
    return rme_target_cover[:, closest_scenario_idx]
end

"""
ReefMod results from 2008 to 2024 follow the same DHW scenarios, the only difference is the
initial cover. Up until 2024, from the 100 scenarios, there are only 20 unique for each
year. From 2025 on we have 100 unique for each year.
"""
function _reefmod_cover(
    reefmod_path::String,
    historic_data::YAXArray{Union{Missing, T},2},
    canonical_gpkg::DataFrame
) where T<:AbstractFloat
    # Get canonical indexes for each historic_data location
    historic_locs = historic_data.locations
    canonical_locs = parse.(Int64, canonical_gpkg.RME_UNIQUE_ID)
    canonical_idx = [findfirst(hist_loc .== canonical_locs) for hist_loc in historic_locs]

    # Return RME coral cover with locations in the same order as in the transect_cube
    reefmod_ds = open_dataset(reefmod_path)
    return reefmod_ds.coral_cover[location=canonical_idx]
end

function _remove_missing_unique_id(ltmp_data)
    return ltmp_data[.!ismissing.(ltmp_data.RME_UNIQUE_ID), :]
end

function _remove_duplicated_locs(ltmp_data)
    ltmp_idx = ltmp_data.RME_UNIQUE_ID
    ltmp_idx_count = countmap(ltmp_idx)

    # Remove duplicated locations
    duplicated_idx = collect(keys(ltmp_idx_count))[values(ltmp_idx_count).!=1]
    return ltmp_data[ltmp_idx.∉[duplicated_idx], :]
end

function _ltmp_cube(ltmp_data_path::String)::YAXArray
    # Load and clean up LTMP manta tow data
    ltmp_gpkg = GDF.read(ltmp_data_path)
    ltmp_gpkg = _remove_missing_unique_id(ltmp_gpkg)
    ltmp_gpkg = _remove_duplicated_locs(ltmp_gpkg)

    ltmp_gpkg_col_names = names(ltmp_gpkg)
    timesteps = sort(ltmp_gpkg_col_names[first.(ltmp_gpkg_col_names, 2) .∈ [["19", "20"]]])
    locations = parse.(Int64, ltmp_gpkg.RME_UNIQUE_ID)
    @assert length(unique(locations)) == length(locations)

    n_timesteps = length(timesteps)
    n_locations = length(locations)
    ltmp_data = zeros(Union{Missing, Float64}, n_timesteps, n_locations)
    for (location_id, location) in enumerate(locations)
        ltmp_data[:,location_id] .= Matrix(ltmp_gpkg[locations .== location, timesteps])[1,:]
    end

    return ADRIA.DataCube(ltmp_data, timesteps=parse.(Int64, timesteps), locations=locations)
end

"""
    _transect_cube(transect_path::String, canonical_gpkg::DataFrame)
    _transect_cube(transect_df::DataFrame)::YAXArray{Float64,2}

# Arguments
- `transect_path` : Path to `LTMP_TRANSECT_MEANS.parquet` file
- `canonical_gpkg` : Path to canonical geopackage
- `transect_df` : LTMP Transect means DataFrame
"""
function _transect_cube(
    transect_path::String, canonical_gpkg::DataFrame
)::YAXArray
    transect_ds = Parquet2.Dataset(transect_path)
    transect_df = DataFrame(transect_ds)
    transect_df.UNIQUE_ID = parse.(Int64, _transect_unique_ids(transect_df, canonical_gpkg))
    return _transect_cube(transect_df)
end
function _transect_cube(transect_df::DataFrame)::YAXArray
    timesteps = Vector{Float64}(sort(unique(transect_df.YEAR_Reefmod)))
    timesteps = Int64.(timesteps[timesteps.% 1 .== 0.0])
    locations = unique(transect_df.UNIQUE_ID)
    n_timesteps = length(timesteps)
    n_locations = length(locations)
    transect_data = zeros(Union{Float64, Missing}, n_timesteps, n_locations)
    for (idx_timestep, timestep) in enumerate(timesteps)
        for (idx_location, location) in enumerate(locations)
            location_filter = transect_df.UNIQUE_ID.==location
            timestep_filter = transect_df.YEAR_Reefmod.==timestep
            timestep_location_filter = timestep_filter.&&location_filter
            cover = Vector{Float64}(transect_df[timestep_location_filter,:Fun_Fun_TOT])
            transect_data[idx_timestep, idx_location] = isempty(cover) ? missing : cover[1]
        end
    end

    return ADRIA.DataCube(transect_data; timesteps=timesteps, locations=locations)
end

function _transect_unique_ids(
    transect_df::DataFrame, canonical_gpkg::DataFrame
)::Vector{String}
    transect_reefs = transect_df.KarloReefName
    canonical_reefs = canonical_gpkg.reef_name
    match_rows = [findfirst(t_reef .== canonical_reefs) for t_reef in transect_reefs]
    return canonical_gpkg[match_rows, "RME_UNIQUE_ID"]
end
