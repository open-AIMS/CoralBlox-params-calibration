# Reads data/clean_growth_intervals.csv and writes data/background_growth_rates.csv: one
# linear model of background (undisturbed) relative per-year coral cover growth, fit on TOTAL
# manta tow cover (not per functional group), as a function of starting cover and survey-gap
# length, with a cluster-bootstrap-by-reef CI. Used by 05_create_survival_rates_csv.jl to
# correct survival rates for background growth. See README.md, Decisions, for the model's
# design history.
#
# Run from the repo root: julia --project=. src/historical_disturbance_mortality_data_gen/02_estimate_background_growth_rates.jl

using CSV, DataFrames, Statistics, Random, LinearAlgebra

DATA_DIR = joinpath(@__DIR__, "data")

# Same fixed seed convention as _per_reef_bootstrap_stats in src/common/perf_metrics.jl, so
# this file's numbers are exactly reproducible from clean_growth_intervals.csv alone.
RNG_SEED = 1
Random.seed!(RNG_SEED)
BOOTSTRAP_REPLICATES = 2000
CI_LEVEL = 0.95

clean_intervals_path = joinpath(DATA_DIR, "clean_growth_intervals.csv")
clean_intervals = CSV.read(clean_intervals_path, DataFrame; stringtype=String, comment="#")
println("Clean interval candidates: ", nrow(clean_intervals))

# Matches the calibration period (START_YEAR in src/common/constants.jl).
MIN_YEAR_START = 2008
if MIN_YEAR_START !== nothing
    n_before = nrow(clean_intervals)
    clean_intervals = clean_intervals[clean_intervals.year_start.>=MIN_YEAR_START, :]
    println("Dropped ", n_before - nrow(clean_intervals), " intervals with year_start < $MIN_YEAR_START")
end

# Avoids division-by-near-zero noise in the relative-growth calculation below.
MIN_COVER_FOR_RELATIVE_GROWTH = 0.005

# relative_growth is an ANNUALIZED (compounding) rate: the constant per-year rate that,
# compounded over gap_years, takes cover_start to cover_end.
growth_observations = DataFrame(
    unique_id=String[], gap_years=Float64[], cover_start=Float64[], relative_growth=Float64[]
)
for interval in eachrow(clean_intervals)
    cover_start = interval.manta_cover_start
    cover_start <= MIN_COVER_FOR_RELATIVE_GROWTH && continue
    cover_end = interval.manta_cover_end
    relative_growth = (max(cover_end, 0.0) / cover_start)^(1 / interval.gap_years) - 1
    push!(growth_observations, (
        string(interval.unique_id), Float64(interval.gap_years), cover_start, relative_growth
    ))
end
println("Growth observations (total-cover intervals): ", nrow(growth_observations))

# Per-row weight `1 / (intervals from that row's reef)`, so each reef contributes equal total
# weight regardless of how many overlapping intervals it happens to supply (reefs with long
# survey runs generate combinatorially more pairwise intervals - see README.md, Decisions).
# Computed once from the unresampled data and carried as a fixed column, so a reef drawn twice
# by the bootstrap below still carries double weight, not weight normalized back to 1.
reef_interval_counts = Dict{String,Int}()
for id in growth_observations.unique_id
    reef_interval_counts[id] = get(reef_interval_counts, id, 0) + 1
end
growth_observations.weight = [1 / reef_interval_counts[id] for id in growth_observations.unique_id]

"""
    design_matrix(rows)::Matrix{Float64}

`[1 cover_start gap_years]` per row - the regressors for `relative_growth ~
intercept + slope_cover_start * cover_start + slope_gap * gap_years`.
"""
design_matrix(rows) = hcat(ones(nrow(rows)), rows.cover_start, rows.gap_years)

"""
    weighted_ols(rows)::Vector{Float64}

WLS fit of `relative_growth ~ intercept + slope_cover_start * cover_start + slope_gap *
gap_years`, weighted by `rows.weight` (see its definition above) - solved via the standard
sqrt-weight scaling trick so it reduces to an ordinary `\\` solve.
"""
function weighted_ols(rows::DataFrame)::Vector{Float64}
    sqrt_w = sqrt.(rows.weight)
    return (sqrt_w .* design_matrix(rows)) \ (sqrt_w .* rows.relative_growth)
end

"""
    cluster_bootstrap_ols_ci(rows; B, ci_level)::NamedTuple

Cluster-by-reef bootstrap CI for the WLS coefficients of `relative_growth ~ intercept +
slope_cover_start * cover_start + slope_gap * gap_years`. Resamples reef identity with
replacement, taking every interval (and its fixed `weight`) from each sampled reef, so a reef
drawn twice contributes double weight to that replicate.
"""
function cluster_bootstrap_ols_ci(
    rows::DataFrame; B::Int=BOOTSTRAP_REPLICATES, ci_level::Float64=CI_LEVEL
)::NamedTuple
    point_estimate = weighted_ols(rows)

    reef_ids = unique(rows.unique_id)
    n_reefs = length(reef_ids)
    rows_by_reef = Dict(id => rows[rows.unique_id.==id, :] for id in reef_ids)

    replicates = Matrix{Float64}(undef, B, 3)
    b = 1
    while b <= B
        sampled_reefs = reef_ids[rand(1:n_reefs, n_reefs)]
        resample = reduce(vcat, (rows_by_reef[id] for id in sampled_reefs))
        coefs = weighted_ols(resample)
        all(isfinite, coefs) || continue
        replicates[b, :] .= coefs
        b += 1
    end

    lo_quantile, hi_quantile = (1 - ci_level) / 2, 1 - (1 - ci_level) / 2
    lo = [quantile(replicates[:, k], lo_quantile) for k in 1:3]
    hi = [quantile(replicates[:, k], hi_quantile) for k in 1:3]
    return (intercept=point_estimate[1], intercept_lo=lo[1], intercept_hi=hi[1],
        slope_cover_start=point_estimate[2], slope_cover_start_lo=lo[2], slope_cover_start_hi=hi[2],
        slope_gap=point_estimate[3], slope_gap_lo=lo[3], slope_gap_hi=hi[3])
end

n_reefs = length(unique(growth_observations.unique_id))
min_reefs_for_bootstrap = 3
@assert n_reefs >= min_reefs_for_bootstrap "Only $n_reefs reef(s) - need >= $min_reefs_for_bootstrap for a bootstrap CI"

fit = cluster_bootstrap_ols_ci(growth_observations)
background_growth_rates = DataFrame(
    intercept=[fit.intercept], intercept_lo=[fit.intercept_lo], intercept_hi=[fit.intercept_hi],
    slope_cover_start=[fit.slope_cover_start], slope_cover_start_lo=[fit.slope_cover_start_lo],
    slope_cover_start_hi=[fit.slope_cover_start_hi],
    slope_gap=[fit.slope_gap], slope_gap_lo=[fit.slope_gap_lo], slope_gap_hi=[fit.slope_gap_hi],
    n_intervals=[nrow(growth_observations)], n_reefs=[n_reefs],
)

out_path = joinpath(DATA_DIR, "background_growth_rates.csv")
open(out_path, "w") do io
    println(io, "# Background (undisturbed) annualized relative per-year coral cover growth, fit on TOTAL")
    println(io, "# manta tow cover: relative_growth = intercept + slope_cover_start*cover_start +")
    println(io, "# slope_gap*gap_years. _lo/_hi = 95% cluster-bootstrap-by-reef CI. See README.md,")
    println(io, "# Decisions, for the model's design history. Used by 05_create_survival_rates_csv.jl.")
    println(io, "# Generated by 02_estimate_background_growth_rates.jl.")
end
CSV.write(out_path, background_growth_rates; append=true, writeheader=true)
println("Written to: ", out_path)
