# Question: does per-year background growth rate vary with survey-gap length, per functional
# group? Motivates gap_years as a predictor in 02_estimate_background_growth_rates.jl's
# regression - see README.md, Decisions.
#
# Run from the repo root: julia --project=. src/historical_disturbance_mortality_data_gen/diagnostics/A_growth_rate_by_gap_length.jl

using DataFrames, Statistics
include(joinpath(@__DIR__, "_shared.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "bootstrap_stats.jl"))

DATA_DIR = joinpath(@__DIR__, "..", "data")
long_format_rows = load_long_format_growth_observations(DATA_DIR)
println("Long-format rows: ", nrow(long_format_rows))

# CI: aggregate each reef to one point estimate (median) per (functional_group, gap) bucket
# first, then bootstrap the median across reefs - respects the fact that intervals from the
# same reef aren't independent draws.
println("\n=== DIAGNOSTIC A: per-year growth rate by gap length (mean [n_intervals, n_reefs]) ===")
for functional_group in FUNCTIONAL_GROUPS
    println("-- $functional_group --")
    group_rows = long_format_rows[long_format_rows.functional_group.==functional_group, :]
    for gap in 1:6
        gap_rows = group_rows[group_rows.gap_years.==gap, :]
        isempty(gap_rows) && continue
        per_reef_median = combine(groupby(gap_rows, :unique_id), :per_year_growth => median => :reef_median)
        if nrow(per_reef_median) >= 3
            ci = bootstrap_median_ci(per_reef_median.reef_median)
            ci_text = " 95% CI=[$(round(ci.lo, digits=5)), $(round(ci.hi, digits=5))]"
        else
            ci_text = " 95% CI=n/a (<3 reefs)"
        end
        println("  gap=$gap: mean=", round(mean(gap_rows.per_year_growth), digits=5),
            " median=", round(median(gap_rows.per_year_growth), digits=5), ci_text,
            " n=", nrow(gap_rows), " reefs=", nrow(per_reef_median))
    end
end
println("\nDone.")
