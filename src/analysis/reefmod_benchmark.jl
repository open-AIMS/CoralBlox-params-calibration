using StatsBase
using Parquet
using YAXArrays, DataFrames, NetCDF
using ADRIA: GDF, DataCube

include("../common/perf_metrics.jl")

valid_historic_data_type() = (:manta_tow, :transect)

"""
    rmse_scores(historic_data_path::String, reefmod_path::String, canonical_path::String, historic_data_type::Symbol)::NamedTuple
    rmse_scores(historic_cover::YAXArray{Union{Missing, T1}, 2},model_cover::YAXArray{T2, 3})::NamedTuple where {T1<:AbstractFloat, T2<:AbstractFloat}

RMSE (Root Mean Square Error) for both ReefMod and LTMP average historic data benchmark.
The error is computed against LTMP historic data (either `:manta_tow` or `"transect`). LTMP
data for locations for which the field RME_UNIQUE_ID was missing and for locations with
duplicated data are discarded. Data for locations for which there is less than 4 years of
data collected are also discarded. For each location this selects the first year for which
there is LTMP data within the analyzed timeframe and looks for the ReefMod scenario with the
closest value for that year; selects only the ReefMod years for which there is LTMP data;
computes the RMSE score for ReefMod data points and LTMP data; and computes the RMSE score
for LTMP the benchmark (average of LTMP data for those years and that location) and LTMP
data.

# Arguments
- `historic_data_path` : Path to historic dataset used for the analysis. When
`historic_data_type` is `:manta_tow`, this should be the path to the manta tow LTMP
geopackage (manta_tow_data_reef_lvl.gpkg); when `historic_data_type` is `:transect` this
should be the path to the transect LTMP dataset (`LTMP_TRANSECT_MEANS.parquet`).
- `reefmod_path` : Path to ReefMod model cover results from 2008 to 2022.
- `canonical_path` : Path to canonical GBR geopackage (`canonical_gbr_2024-04-23.gpkg`).
- `historic_data_type` : Type of historic data used for the analysis. Valid options are
$(valid_historic_data_type())
- `historic_data` : YAXArray with historic data with dimensions `(timesteps ⋅ locations)`.
- `model_cover` : YAXArray with ReefMod model cover with dimensions
`(timestep ⋅ location ⋅ scenario)`.

# Returns
NamedTuple with three elements: rme_scores (vector of RMSE scores for ReefMod model results
and historic data for each location), benchmark_scores (vector of RMSE scores for
benchmark, historic average, and historic data), n_observations (vector with number of
observations for each location)

# Example
```
reefmod_path = "datasets/rme_data/rme_cover_20_reps_2008_2022.nc"
canonical_path = "datasets/spatial_data/canonical_gbr_2024-04-23.gpkg"

transect_data_path = "datasets/ltmp_data/LTMP_TRANSECT_MEANS.parquet"
transect_data_type = :transect
model_scores_transect, benchmark_scores_transect, n_observations_transect = rmse_scores(
    transect_data_path, reefmod_path, canonical_path, transect_data_type
)

manta_tow_data_path = "datasets/ltmp_data/manta_tow_data_reef_lvl.gpkg"
manta_tow_data_type = :manta_tow
model_scores_manta_tow, benchmark_scores_manta_tow, n_observations_manta_tow = rmse_scores(
    manta_tow_data_path, reefmod_path, canonical_path, manta_tow_data_type
)
```
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
    historic_function::Function = eval(Meta.parse("_$(historic_data_type)_cube"))
    historic_cover::YAXArray = historic_function(historic_data_path, canonical_gpkg)

    # The cover used here comes from ReefMod model, but that can change in the future
    model_cover::YAXArray = _reefmod_cover(reefmod_path, historic_cover, canonical_gpkg)
    return rmse_scores(historic_cover, model_cover)
end
function rmse_scores(
    historic_cover::YAXArray{Union{Missing, T1}, 2},
    model_cover::YAXArray{T2, 3}
)::NamedTuple where {T1<:AbstractFloat, T2<:AbstractFloat}
    target_timesteps::Vector{Int64} = _intersect_timerange(historic_cover, model_cover)
    historic_cover = _filter_low_data_locations(historic_cover[timesteps=At(target_timesteps)])

    n_locations = size(historic_cover, :locations)
    model_scores::Vector{Float64} = zeros(Float64, n_locations)
    benchmark_scores::Vector{Float64} = zeros(Float64, n_locations)
    n_observations::Vector{Int64} = zeros(Float64, n_locations)
    for (loc_id, historic_cover_row) in enumerate(eachcol(historic_cover))
        # Identify historic years with data
        years_with_data::BitVector = .!ismissing.(collect(historic_cover_row))
        historic_target_cover = @view(historic_cover_row[years_with_data])
        target_years = @view(target_timesteps[years_with_data])

        # Fill scores
        model_scores[loc_id] = _model_scores(model_cover, historic_target_cover, target_years, loc_id)
        benchmark_scores[loc_id] = _benchmark_scores(historic_target_cover)
        n_observations[loc_id] = length(historic_target_cover)
    end
    return (
        model_scores=model_scores,
        benchmark_scores=benchmark_scores,
        n_observations=n_observations
    )
end

function _intersect_timerange(
    historic_cover::YAXArray{Union{Missing, T1}, 2},
    model_cover::YAXArray{T2, 3}
)::Vector{Int64} where {T1<:AbstractFloat, T2<:AbstractFloat}
    start_year = max(minimum(historic_cover.timesteps), minimum(model_cover.timestep))
    end_year = min(maximum(historic_cover.timesteps), maximum(model_cover.timestep))
    return start_year:end_year
end

function _filter_low_data_locations(
    historic_data::YAXArray{Union{Missing, T}, 2}; min_n_data::Int64=4
)::YAXArray{Union{Missing, T}, 2} where {T<:AbstractFloat}
    has_data = .!ismissing.(historic_data)
    target_locs = dropdims(sum(has_data, dims=:timesteps) .>= min_n_data, dims=:timesteps)
    return historic_data[locations=target_locs]
end

function _benchmark_scores(
    historic_target_cover::YAXArray{Union{Missing, T}, 1}
)::T where {T<:AbstractFloat}
    benchmark_cover = fill(mean(historic_target_cover), length(historic_target_cover))
    return rmse(benchmark_cover, historic_target_cover.data)
end

function _model_scores(
    model_cover::YAXArray{T1, 3},
    historic_target_cover::YAXArray{Union{Missing, T2}},
    target_years::Union{SubArray{T3}, Vector{T3}},
    loc_id::T3
)::Float64 where {T1<:AbstractFloat, T2<:AbstractFloat, T3<:Integer}
    model_target_cover = @view(model_cover[timestep=At(target_years), location=loc_id])
    model_target_cover = _model_closest_scenario(model_target_cover, historic_target_cover)
    rmse(model_target_cover, historic_target_cover.data)
end

"""
Find scenario in model_target_cover whose cover at the first year is closest to the
historic_target_cover first year cover.
"""
function _model_closest_scenario(
    model_target_cover::YAXArray{T1, 2},
    historic_target_cover::YAXArray{Union{Missing, T2}, 1}
)::YAXArray{T1, 1} where {T1<:AbstractFloat, T2<:AbstractFloat}
    historic_first_year_cover::T2 = first(historic_target_cover)
    closest_scenario_idx::Int64 = argmin(abs.(historic_first_year_cover .- model_target_cover[1, :]))
    return model_target_cover[:, closest_scenario_idx]
end

"""
    _reefmod_cover(reefmod_path::String,historic_data::YAXArray{Union{Missing, T},2},canonical_gpkg::DataFrame) where {T<:AbstractFloat}

ReefMod results from 2008 to 2024 follow the same DHW scenarios, the only difference is the
initial cover. Up until 2024, from the 100 scenarios, there are only 20 unique for each
year. From 2025 on we have 100 unique for each year.
"""
function _reefmod_cover(
    reefmod_path::String,
    historic_data::YAXArray{Union{Missing, T},2},
    canonical_gpkg::DataFrame
) where {T<:AbstractFloat}
    # Get canonical row indexes for each historic_data location
    historic_locs = historic_data.locations
    canonical_locs = parse.(Int64, canonical_gpkg.RME_UNIQUE_ID)
    canonical_row_idx = [findfirst(hist_loc .== canonical_locs) for hist_loc in historic_locs]

    # Return ReefMod coral cover with locations in the same order as in the historic_data
    reefmod_ds = open_dataset(reefmod_path)
    return reefmod_ds.coral_cover[location=canonical_row_idx]
end

function _remove_missing_unique_id(manta_tow_data::DataFrame)::DataFrame
    return manta_tow_data[.!ismissing.(manta_tow_data.RME_UNIQUE_ID), :]
end

function _remove_duplicated_locs(manta_tow_data::DataFrame)::DataFrame
    manta_tow_idx = manta_tow_data.RME_UNIQUE_ID
    manta_tow_idx_count = countmap(manta_tow_idx)

    # Remove duplicated locations
    duplicated_idx = collect(keys(manta_tow_idx_count))[values(manta_tow_idx_count).!=1]
    return manta_tow_data[manta_tow_idx.∉[duplicated_idx], :]
end

"""
    _manta_tow_cube(manta_tow_data_path::String, _canonical_gpkg::DataFrame)::YAXArray{Union{Missing, Float64}, 2}

# Arguments
- `manta_tow_data_path` : Path to LTMP manta tow dataset.
- `_canonical_gpkg` : Canonical geopackage.
"""
function _manta_tow_cube(
    manta_tow_data_path::String, _canonical_gpkg::DataFrame
)::YAXArray{Union{Missing, Float64}, 2}
    # Load and clean up LTMP manta tow data
    manta_tow_gpkg = GDF.read(manta_tow_data_path)
    manta_tow_gpkg = _remove_missing_unique_id(manta_tow_gpkg)
    manta_tow_gpkg = _remove_duplicated_locs(manta_tow_gpkg)

    manta_tow_col_names = names(manta_tow_gpkg)
    timesteps = sort(manta_tow_col_names[first.(manta_tow_col_names, 2) .∈ [["19", "20"]]])
    locations = parse.(Int64, manta_tow_gpkg.RME_UNIQUE_ID)
    @assert length(unique(locations)) == length(locations)

    n_timesteps = length(timesteps)
    n_locations = length(locations)
    manta_tow_data = zeros(Union{Missing, Float64}, n_timesteps, n_locations)
    for (location_id, location) in enumerate(locations)
        manta_tow_data[:,location_id] .= Matrix(manta_tow_gpkg[locations .== location, timesteps])[1,:]
    end

    return DataCube(manta_tow_data, timesteps=parse.(Int64, timesteps), locations=locations)
end

"""
    _transect_cube(transect_path::String, canonical_gpkg::DataFrame)::YAXArray{Union{Missing, Float64}, 2}
    _transect_cube(transect_df::DataFrame)::YAXArray{Union{Missing, Float64}, 2}

# Arguments
- `transect_path` : Path to `LTMP_TRANSECT_MEANS.parquet` file
- `canonical_gpkg` : Path to canonical geopackage
- `transect_df` : LTMP Transect means DataFrame
"""
function _transect_cube(
    transect_path::String, canonical_gpkg::DataFrame
)::YAXArray{Union{Missing, Float64}, 2}
    transect_ds = Parquet.read_parquet(transect_path)
    transect_df = DataFrame(transect_ds)
    transect_df.UNIQUE_ID = parse.(Int64, _transect_unique_ids(transect_df, canonical_gpkg))
    return _transect_cube(transect_df)
end
function _transect_cube(
    transect_df::DataFrame
)::YAXArray{Union{Missing, Float64}, 2}
    # Select only whole years as timesteps (remove 2000.5, 2001.5, etc)
    timesteps = Vector{Float64}(sort(unique(transect_df.YEAR_Reefmod)))
    timesteps = Int64.(timesteps[timesteps.% 1 .== 0.0])
    n_timesteps::Int64 = length(timesteps)

    locations = unique(transect_df.UNIQUE_ID)
    n_locations::Int64 = length(locations)
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

    return DataCube(transect_data; timesteps=timesteps, locations=locations)
end

"""
    For each row in `transect_df`, find the correspondent `canonical_gpkg`'s "RME_UNIQUE_ID"
that matches that row's reef.
"""
function _transect_unique_ids(
    transect_df::DataFrame, canonical_gpkg::DataFrame
)::Vector{String}
    transect_reefs = transect_df.KarloReefName
    canonical_reefs = canonical_gpkg.reef_name
    reef_row_indexes = [findfirst(t_reef .== canonical_reefs) for t_reef in transect_reefs]
    return canonical_gpkg[reef_row_indexes, "RME_UNIQUE_ID"]
end
