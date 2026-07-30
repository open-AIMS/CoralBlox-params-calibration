# Question: are gap length and starting cover themselves correlated? Diagnostic A attributes
# curvature to gap length and Diagnostic B/C attribute (weak, or nonlinear) effects to
# starting cover, but those were tested independently - if gap length and starting cover are
# correlated (e.g. long-gap intervals disproportionately come from low-cover reefs), part of
# the "gap effect" could actually be a cover effect in disguise, or vice versa.
#
# Run from the repo root: julia --project=. src/historical_disturbance_mortality_data_gen/diagnostics/D_gap_cover_confound.jl

using DataFrames, Statistics
include(joinpath(@__DIR__, "_shared.jl"))

DATA_DIR = joinpath(@__DIR__, "..", "data")
long_format_rows = load_long_format_growth_observations(DATA_DIR)
println("Long-format rows: ", nrow(long_format_rows))

println("\n=== DIAGNOSTIC D: gap_years vs cover_start correlation (confound check) ===")
for functional_group in FUNCTIONAL_GROUPS
    group_rows = long_format_rows[long_format_rows.functional_group.==functional_group, :]
    correlation = cor(group_rows.gap_years, group_rows.cover_start)
    println("  $functional_group: r=", round(correlation, digits=3), " n=", nrow(group_rows))
end
println("\nDone.")
