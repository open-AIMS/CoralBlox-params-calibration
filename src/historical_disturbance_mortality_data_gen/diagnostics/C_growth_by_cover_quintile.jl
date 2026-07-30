# Question: does binning by starting-cover quintile reveal a nonlinear density-dependence
# effect that Diagnostic B's linear correlation could hide (e.g. growth staying flat until
# cover approaches a ceiling, then decelerating only in the top bin(s))?
#
# Run from the repo root: julia --project=. src/historical_disturbance_mortality_data_gen/diagnostics/C_growth_by_cover_quintile.jl

using DataFrames, Statistics
include(joinpath(@__DIR__, "_shared.jl"))

DATA_DIR = joinpath(@__DIR__, "..", "data")
long_format_rows = load_long_format_growth_observations(DATA_DIR)
println("Long-format rows: ", nrow(long_format_rows))

println("\n=== DIAGNOSTIC C: mean per-year growth rate by starting-cover quintile (absolute) ===")
n_quintiles = 5
for functional_group in FUNCTIONAL_GROUPS
    println("-- $functional_group --")
    group_rows = long_format_rows[long_format_rows.functional_group.==functional_group, :]
    quintile_edges = quantile(group_rows.cover_start, range(0, 1; length=n_quintiles + 1))
    for quintile in 1:n_quintiles
        lower_edge, upper_edge = quintile_edges[quintile], quintile_edges[quintile + 1]
        in_quintile = quintile == n_quintiles ?
                      (group_rows.cover_start .>= lower_edge) .& (group_rows.cover_start .<= upper_edge) :
                      (group_rows.cover_start .>= lower_edge) .& (group_rows.cover_start .< upper_edge)
        quintile_rows = group_rows[in_quintile, :]
        isempty(quintile_rows) && continue
        println("  quintile $quintile  cover_start ∈ [", round(lower_edge, digits=4), ", ", round(upper_edge, digits=4), ")",
            "  mean_growth=", round(mean(quintile_rows.per_year_growth), digits=5),
            "  n=", nrow(quintile_rows))
    end
end

# Same binning, on the RELATIVE growth measure production actually uses.
println("\n=== DIAGNOSTIC C (relative): mean per-year relative growth rate by starting-cover quintile ===")
for functional_group in FUNCTIONAL_GROUPS
    println("-- $functional_group --")
    group_rows = long_format_rows[long_format_rows.functional_group.==functional_group, :]
    group_rows = group_rows[.!ismissing.(group_rows.per_year_relative_growth), :]
    relative_growth = Float64.(group_rows.per_year_relative_growth)
    quintile_edges = quantile(group_rows.cover_start, range(0, 1; length=n_quintiles + 1))
    for quintile in 1:n_quintiles
        lower_edge, upper_edge = quintile_edges[quintile], quintile_edges[quintile + 1]
        in_quintile = quintile == n_quintiles ?
                      (group_rows.cover_start .>= lower_edge) .& (group_rows.cover_start .<= upper_edge) :
                      (group_rows.cover_start .>= lower_edge) .& (group_rows.cover_start .< upper_edge)
        isempty(relative_growth[in_quintile]) && continue
        println("  quintile $quintile  cover_start ∈ [", round(lower_edge, digits=4), ", ", round(upper_edge, digits=4), ")",
            "  mean_relative_growth=", round(mean(relative_growth[in_quintile]), digits=4),
            "  median_relative_growth=", round(median(relative_growth[in_quintile]), digits=4),
            "  n=", count(in_quintile))
    end
end
println("\nDone.")
