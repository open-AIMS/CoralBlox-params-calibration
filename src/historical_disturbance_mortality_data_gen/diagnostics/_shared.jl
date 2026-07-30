# Shared helper for the scripts in this folder - not a diagnostic itself. Reshapes
# data/clean_growth_intervals.csv into one row per (interval, functional_group), with
# per-year growth in both absolute (per_year_growth) and relative (per_year_relative_growth)
# terms - see README.md, Decisions, for why the relative measure is what production uses.

using CSV, DataFrames

# isdefined-guarded so this file can be `include`d more than once in the same session (e.g.
# running two diagnostic scripts back to back) without erroring on a duplicate `const`.
isdefined(Main, :FUNCTIONAL_GROUPS) || (const FUNCTIONAL_GROUPS = [
    "Tabular_Acropora", "Corymbose_Acropora", "Corymbose_non_Acropora",
    "Small_Massives", "Large_Massives"
])

# Avoids division-by-near-zero noise in the relative growth measure.
isdefined(Main, :MIN_COVER_FOR_RELATIVE_GROWTH) || (const MIN_COVER_FOR_RELATIVE_GROWTH = 0.005)

function load_long_format_growth_observations(data_dir::AbstractString)::DataFrame
    clean_intervals_path = joinpath(data_dir, "clean_growth_intervals.csv")
    clean_intervals = CSV.read(clean_intervals_path, DataFrame; stringtype=String, comment="#")

    long_format_rows = DataFrame(
        unique_id=String[], reef_name=String[], year_start=Int[], year_end=Int[], gap_years=Int[],
        functional_group=String[], cover_start=Float64[], cover_end=Float64[],
        per_year_growth=Float64[], per_year_relative_growth=Union{Missing,Float64}[]
    )
    for interval in eachrow(clean_intervals)
        for functional_group in FUNCTIONAL_GROUPS
            cover_start = interval[Symbol("cover_start_$functional_group")]
            cover_end = interval[Symbol("cover_end_$functional_group")]
            relative_growth = cover_start > MIN_COVER_FOR_RELATIVE_GROWTH ?
                               (cover_end - cover_start) / (cover_start * interval.gap_years) :
                               missing
            push!(long_format_rows, (
                string(interval.unique_id), interval.reef_name, interval.year_start, interval.year_end,
                interval.gap_years, functional_group, cover_start, cover_end,
                (cover_end - cover_start) / interval.gap_years, relative_growth
            ))
        end
    end
    return long_format_rows
end
