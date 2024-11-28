using StatsBase
using YAXArrays, DataFrames, NetCDF
using ADRIA: GDF

include("../common/perf_metrics.jl")

function rmse_scores(ltmp_data_path, rme_benchmark_path, canonical_path)
    # Load LTMP and RME data
    start_year, end_year = 2008, 2022
    _ltmp_data::DataFrame = ltmp_data(
        ltmp_data_path;
        start_year=start_year, end_year=end_year
    )
    _rme_cover::YAXArray = rme_cover(
        Vector{String}(_ltmp_data.RME_UNIQUE_ID), rme_benchmark_path, canonical_path;
        start_year=start_year, end_year=end_year
    )

    return rmse_scores(_ltmp_data, _rme_cover)
end

function rmse_scores(ltmp_data, rme_cover)
    n_locations::Int64 = size(ltmp_data, 1)
    rme_scores::Vector{Float64} = zeros(Float64, n_locations)
    benchmark_scores::Vector{Float64} = zeros(Float64, n_locations)
    n_observations::Vector{Float64} = zeros(Float64, n_locations)

    start_year, end_year = extrema(rme_cover.timestep)
    for (row_idx, ltmp_row) in enumerate(eachrow(ltmp_data[:, string.(start_year:end_year)]))
        # Identify LTMP years with data
        years_with_data = .!ismissing.(collect(ltmp_row))
        ltmp_row_with_data = ltmp_row[years_with_data]
        n_observations[row_idx] = length(ltmp_row_with_data)

        first_year = findmin(names(ltmp_row_with_data))[1]
        first_year_ltmp_data = ltmp_row_with_data[first_year]

        # Find scenario with closest cover for that year
        rme_cover_scenarios = rme_cover[timestep=years_with_data, location=row_idx]
        rme_unique_scens = hcat(unique(eachcol(rme_cover_scenarios.data))...)
        target_scen_idx = argmin(abs.(first_year_ltmp_data .- rme_unique_scens[1, :]))

        rme_y = rme_unique_scens[:, target_scen_idx]
        ltmp_y = collect(ltmp_row_with_data)
        benchmark_y = fill(mean(ltmp_y), length(ltmp_y))

        # Compute scores
        rme_scores[row_idx] = rmse(rme_y, ltmp_y)
        benchmark_scores[row_idx] = rmse(benchmark_y, ltmp_y)
    end

    (
        rme_scores=rme_scores,
        benchmark_scores=benchmark_scores,
        n_observations=n_observations
    )
end

function ltmp_data(
    ltmp_data_path::String;
    min_n_data::Int64=4, start_year::Int64=2008, end_year::Int64=2022
)::DataFrame
    # Load and clean up LTMP manta tow data
    _ltmp_data = GDF.read(ltmp_data_path)
    _ltmp_data = _remove_missing_unique_id(_ltmp_data)
    _ltmp_data = _remove_duplicated_locs(_ltmp_data)

    n_data_per_loc = dropdims(sum(Matrix(.!ismissing.(_ltmp_data[:, string.(start_year:end_year)])), dims=2), dims=2)
    return _ltmp_data[n_data_per_loc.>=min_n_data, :]
end

"""
ReefMod results from 2008 to 2024 follow the same DHW scenarios, the only difference is the
initial cover. Up until 2024, from the 100 scenarios, there are only 20 unique for each
year. From 2025 on we have 100 unique for each year.
"""
function rme_cover(
    ltmp_idx::Vector{String},
    rme_benchmark_path::String,
    canonical_path::String;
    start_year::Int64=2008,
    end_year::Int64=2022
)::YAXArray
    rme_ds = open_dataset(rme_benchmark_path)

    # Select RME cover only for locations and timesteps for which we have LTMP data
    canonical_idx = GDF.read(canonical_path).UNIQUE_ID
    locations_mask = canonical_idx .∈ [ltmp_idx]
    rme_cover = rme_ds.coral_cover_per_taxa[
        location=locations_mask,
        timestep=At(collect(start_year:end_year))
    ]

    # Sum across all functional groups
    return dropdims(sum(rme_cover, dims=:group), dims=:group)
end

function _remove_missing_unique_id(ltmp_data)
    return ltmp_data[.!ismissing.(ltmp_data.RME_UNIQUE_ID), :]
end

function _remove_duplicated_locs(ltmp_data)
    ltmp_idx = ltmp_data.RME_UNIQUE_ID
    ltmp_idx_count = countmap(ltmp_idx)

    # Duplicated idx "11184103104", "11222101104", "16013101104"
    duplicated_idx = collect(keys(ltmp_idx_count))[values(ltmp_idx_count).!=1]
    return ltmp_data[ltmp_idx.∉[duplicated_idx], :]
end
