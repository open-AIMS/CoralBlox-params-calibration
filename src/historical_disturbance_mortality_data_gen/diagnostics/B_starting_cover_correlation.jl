# Question: is starting cover linearly correlated with per-year growth (a simple density-
# dependence check)? A weak result here doesn't rule out a nonlinear effect - see
# C_growth_by_cover_quintile.jl for that follow-up.
#
# Run from the repo root: julia --project=. src/historical_disturbance_mortality_data_gen/diagnostics/B_starting_cover_correlation.jl

using DataFrames, Statistics
include(joinpath(@__DIR__, "_shared.jl"))

DATA_DIR = joinpath(@__DIR__, "..", "data")
long_format_rows = load_long_format_growth_observations(DATA_DIR)
println("Long-format rows: ", nrow(long_format_rows))

println("\n=== DIAGNOSTIC B: cover_start vs per_year_growth (absolute) correlation ===")
for functional_group in FUNCTIONAL_GROUPS
    group_rows = long_format_rows[long_format_rows.functional_group.==functional_group, :]
    correlation = cor(group_rows.cover_start, group_rows.per_year_growth)
    println("  $functional_group: r=", round(correlation, digits=3), " n=", nrow(group_rows),
        " mean_start_cover=", round(mean(group_rows.cover_start), digits=4))
end

# Same check for the RELATIVE growth measure production actually uses - see README.md,
# Decisions.
println("\n=== DIAGNOSTIC B (relative): cover_start vs per_year_relative_growth correlation ===")
for functional_group in FUNCTIONAL_GROUPS
    group_rows = long_format_rows[long_format_rows.functional_group.==functional_group, :]
    group_rows = group_rows[.!ismissing.(group_rows.per_year_relative_growth), :]
    correlation = cor(group_rows.cover_start, Float64.(group_rows.per_year_relative_growth))
    println("  $functional_group: r=", round(correlation, digits=3), " n=", nrow(group_rows),
        " mean_start_cover=", round(mean(group_rows.cover_start), digits=4))
end
println("\nDone.")
