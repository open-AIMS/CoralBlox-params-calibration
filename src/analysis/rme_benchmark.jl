using StatsBase
using YAXArrays, DataFrames, NetCDF
using ADRIA: GDF

include("../common/perf_metrics.jl")

# TODO write doctring explaining the analysis
"""
    rmse_scores(ltmp_data_path, rme_benchmark_path, canonical_path)
"""
function rmse_scores(
    ltmp_data_path::String,
    rme_benchmark_path::String,
    canonical_path::String
)::NamedTuple
    # Load LTMP and RME data
    start_year, end_year = 2008, 2022
    ltmp_data::DataFrame = _ltmp_data(
        ltmp_data_path; start_year=start_year, end_year=end_year
    )
    ltmp_idx::Vector{String} = ltmp_data.RME_UNIQUE_ID
    rme_cover::YAXArray = _rme_cover(ltmp_idx, rme_benchmark_path, canonical_path)

    ltmp_cover = ltmp_data[:, string.(start_year:end_year)]
    return rmse_scores(ltmp_cover, rme_cover)
end

function rmse_scores(ltmp_cover::DataFrame, rme_cover::YAXArray)::NamedTuple
    n_locations::Int64 = size(ltmp_cover, 1)
    rme_scores::Vector{Float64} = zeros(Float64, n_locations)
    benchmark_scores::Vector{Float64} = zeros(Float64, n_locations)
    n_observations::Vector{Int64} = zeros(Float64, n_locations)
    for (loc_idx, ltmp_loc_cover) in enumerate(eachrow(ltmp_cover))
        # Identify LTMP years with data
        years_with_data::BitVector = .!ismissing.(collect(ltmp_loc_cover))
        ltmp_target_cover = ltmp_loc_cover[years_with_data]
        target_years = parse.(Int64, names(ltmp_target_cover))
        ltmp_raw_cover = collect(ltmp_target_cover)

        # RME scores
        rme_target_cover = rme_cover[timestep=At(target_years), location=loc_idx]
        rme_target_cover = _rme_year_slice(rme_target_cover, ltmp_target_cover)
        rme_scores[loc_idx] = rmse(rme_target_cover, ltmp_raw_cover)

        # Benchmark scores
        n_observations[loc_idx] = length(ltmp_raw_cover)
        benchmark_cover = fill(mean(ltmp_target_cover), n_observations[loc_idx])
        benchmark_scores[loc_idx] = rmse(benchmark_cover, ltmp_raw_cover)
    end
    return (
        rme_scores=rme_scores,
        benchmark_scores=benchmark_scores,
        n_observations=n_observations
    )
end

function _rme_year_slice(rme_loc_cover::YAXArray, ltmp_loc_cover::DataFrameRow)::YAXArray
    first_year_with_ltmp_data = findmin(names(ltmp_loc_cover))[1]
    ltmp_first_year_cover = ltmp_loc_cover[first_year_with_ltmp_data]
    closest_scenario_idx = argmin(abs.(ltmp_first_year_cover .- rme_loc_cover[1, :]))
    return rme_loc_cover[:, closest_scenario_idx]
end

function _ltmp_data(
    ltmp_data_path::String;
    min_n_data::Int64=4, start_year::Int64=2008, end_year::Int64=2022
)::DataFrame
    # Load and clean up LTMP manta tow data
    _ltmp_data = GDF.read(ltmp_data_path)
    _ltmp_data = _remove_missing_unique_id(_ltmp_data)
    _ltmp_data = _remove_duplicated_locs(_ltmp_data)

    # Select LTMP DataFrame cols with cover data from `start_year` to `end_year` to
    # compute number of datapoints for each location
    ltmp_data_years = Matrix(.!ismissing.(_ltmp_data[:, string.(start_year:end_year)]))
    n_data_per_loc = dropdims(sum(ltmp_data_years, dims=2), dims=2)

    return _ltmp_data[n_data_per_loc.>=min_n_data, :]
end

"""
ReefMod results from 2008 to 2024 follow the same DHW scenarios, the only difference is the
initial cover. Up until 2024, from the 100 scenarios, there are only 20 unique for each
year. From 2025 on we have 100 unique for each year.
"""
function _rme_cover(
    ltmp_idx::Vector{String},
    rme_benchmark_path::String,
    canonical_path::String;
)::YAXArray
    rme_ds = open_dataset(rme_benchmark_path)

    # Select RME cover only for locations and timesteps for which we have LTMP data
    canonical_idx = GDF.read(canonical_path).UNIQUE_ID
    locations_mask = canonical_idx .∈ [ltmp_idx]
    return rme_ds.coral_cover[location=locations_mask]
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
