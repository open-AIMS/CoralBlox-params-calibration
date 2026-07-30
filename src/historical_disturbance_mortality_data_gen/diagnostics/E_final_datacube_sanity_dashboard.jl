# Gut-feel sanity dashboard for the final datacube, not an exhaustive validation suite:
#   1. Datacube mortality by functional group (real disturbance events only) - biological
#      plausibility check.
#   2. Intermediate survival_rate distribution, with the [0.05, 1.0] clamp bounds marked.
#   3. CSV-implied vs datacube mortality - pipeline wiring check (should sit on the 1:1 line).
#
# Run from the repo root: julia --project=. src/historical_disturbance_mortality_data_gen/diagnostics/E_final_datacube_sanity_dashboard.jl

using CSV, DataFrames, YAXArrays, NetCDF, CairoMakie, ColorSchemes, TOML
import GeoDataFrames as GDF

include(joinpath(@__DIR__, "_shared.jl"))

DATA_DIR = joinpath(@__DIR__, "..", "data")
REPO_ROOT = joinpath(@__DIR__, "..", "..", "..")
OUTPUT_DIR = joinpath(@__DIR__, "output")
mkpath(OUTPUT_DIR)

# Positionally aligned with FUNCTIONAL_GROUPS (Title_Case, CSV columns) - lower_case, matches
# the NetCDF `species` axis.
SPECIES_NAMES = [
    "tabular_Acropora", "corymbose_Acropora", "corymbose_non_Acropora",
    "small_massives", "large_massives"
]
GROUP_TO_SPECIES = Dict(FUNCTIONAL_GROUPS .=> SPECIES_NAMES)
GROUP_COLORS = Dict(FUNCTIONAL_GROUPS .=> ColorSchemes.seaborn_colorblind6[1:5])

config = TOML.parsefile(joinpath(REPO_ROOT, "config.toml"))
canonical_path = config["calibration"]["geospatial"]["canonical_path"]
canonical_gpkg = GDF.read(canonical_path)
reef_id_to_gbrmpa_id = Dict(canonical_gpkg.RME_UNIQUE_ID .=> canonical_gpkg.RME_GBRMPA_ID)

survival_rates_df = CSV.read(
    joinpath(DATA_DIR, "disturbance_survival_rates.csv"), DataFrame;
    stringtype=String, comment="#"
)

datacube_path = joinpath(DATA_DIR, "historical_disturbance_mortality_rates", "historical_disturbance_mortality_rates.nc")
datacube = open_dataset(datacube_path).disturbance_mortality_scens

# One row per (disturbance event, functional group). dist_year clamped to 2008 to match
# 06_build_dist_mortality_datacube.jl's handling of pre-2008 disturbances.
comparison_rows = DataFrame(
    reef_name=String[], dist_year=Int[], functional_group=String[],
    survival_rate=Float64[], csv_mortality=Float64[], datacube_mortality=Float64[]
)
for disturbance_event in eachrow(survival_rates_df)
    gbrmpa_id = reef_id_to_gbrmpa_id[string(disturbance_event.reef_id)]
    datacube_year = max(disturbance_event.dist_year, 2008)
    for functional_group in FUNCTIONAL_GROUPS
        survival_rate = disturbance_event["survival_rate_$functional_group"]
        csv_mortality = clamp(1 - survival_rate, 0.0, 0.95)
        datacube_mortality = datacube[
            locs=At(gbrmpa_id), timesteps=At(datacube_year),
            species=At(GROUP_TO_SPECIES[functional_group]), scenarios=At(1)
        ][1]
        push!(comparison_rows, (
            disturbance_event.reef_name, disturbance_event.dist_year, functional_group,
            survival_rate, csv_mortality, datacube_mortality
        ))
    end
end
println("Disturbance-event x functional-group rows compared: ", nrow(comparison_rows))

println("\n=== DIAGNOSTIC E: survival_rate clamp hit-rate by functional group ===")
for functional_group in FUNCTIONAL_GROUPS
    group_rows = comparison_rows[comparison_rows.functional_group.==functional_group, :]
    n_total = nrow(group_rows)
    n_floor_clamped = sum(group_rows.survival_rate .<= 0.0500001)
    n_ceiling_clamped = sum(group_rows.survival_rate .>= 0.999999)
    println(
        "  $functional_group: n=$n_total  floor(0.05, near-total mortality)=$n_floor_clamped ",
        "(", round(100 * n_floor_clamped / n_total, digits=1), "%)  ",
        "ceiling(1.0, zero mortality)=$n_ceiling_clamped ",
        "(", round(100 * n_ceiling_clamped / n_total, digits=1), "%)"
    )
end

println("\n=== DIAGNOSTIC E: CSV-implied vs datacube mortality wiring check ===")
max_abs_diff = maximum(abs.(comparison_rows.csv_mortality .- comparison_rows.datacube_mortality))
println("  max |csv_mortality - datacube_mortality| across all rows: ", max_abs_diff)

# ----- Figure -----------------------------------------------------------------

fig = Figure(size=(1500, 500))

ax1 = Axis(
    fig[1, 1], title="Datacube mortality by functional group\n(real disturbance events only)",
    xlabel="Functional group", ylabel="Mortality rate (as stored in NetCDF)",
    xticks=(1:5, replace.(FUNCTIONAL_GROUPS, "_" => "\n")), xticklabelsize=11
)
for (group_index, functional_group) in enumerate(FUNCTIONAL_GROUPS)
    group_rows = comparison_rows[comparison_rows.functional_group.==functional_group, :]
    boxplot!(
        ax1, fill(group_index, nrow(group_rows)), group_rows.datacube_mortality;
        color=GROUP_COLORS[functional_group], show_outliers=true
    )
end

ax2 = Axis(
    fig[1, 2], title="Intermediate survival_rate distribution\n(clamp bounds marked)",
    xlabel="Functional group", ylabel="survival_rate (pre mortality-clamp)",
    xticks=(1:5, replace.(FUNCTIONAL_GROUPS, "_" => "\n")), xticklabelsize=11
)
for (group_index, functional_group) in enumerate(FUNCTIONAL_GROUPS)
    group_rows = comparison_rows[comparison_rows.functional_group.==functional_group, :]
    boxplot!(
        ax2, fill(group_index, nrow(group_rows)), group_rows.survival_rate;
        color=GROUP_COLORS[functional_group], show_outliers=true
    )
end
hlines!(ax2, [0.05, 1.0]; color=:black, linestyle=:dash, linewidth=1)
text!(ax2, 0.55, 0.05, text="floor", fontsize=10, align=(:left, :top))
text!(ax2, 0.55, 1.0, text="ceiling", fontsize=10, align=(:left, :bottom))

ax3 = Axis(
    fig[1, 3], title="Wiring check: CSV-implied vs datacube mortality\n(should sit on 1:1 line)",
    xlabel="Mortality implied by disturbance_survival_rates.csv", ylabel="Mortality in historical_disturbance_mortality_rates.nc"
)
ablines!(ax3, 0, 1; color=:grey, linestyle=:dash)
for functional_group in FUNCTIONAL_GROUPS
    group_rows = comparison_rows[comparison_rows.functional_group.==functional_group, :]
    scatter!(
        ax3, group_rows.csv_mortality, group_rows.datacube_mortality;
        color=GROUP_COLORS[functional_group], label=functional_group, markersize=6
    )
end
axislegend(ax3; position=:lt, labelsize=9, framevisible=false)

output_path = joinpath(OUTPUT_DIR, "E_final_datacube_sanity_dashboard.png")
save(output_path, fig)
println("\nSaved dashboard to: ", output_path)
println("Done.")
