# 01_create_disturbance_years_csv.jl is deliberately NOT in this list - see its own header
# comment for why. Run it standalone if needed.
script_names = [
    "base.jl",
    "02_estimate_background_growth_rates.jl",
    "03_a_create_cover_manta_csv.jl",
    "03_b_create_cover_transect_csv.jl",
    "04_fill_missing_transect_disturbances.jl",
    "05_create_survival_rates_csv.jl",
    "06_build_dist_mortality_datacube.jl",
    "07_create_unmapped_disturbance_years.jl"
]

for name in script_names
    @info "Running $(name)"
    include(name)
end
