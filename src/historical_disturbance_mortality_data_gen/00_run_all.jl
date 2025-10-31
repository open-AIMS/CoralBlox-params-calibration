include("disturbance_data_gen/base.jl")

include("disturbance_data_gen/01_create_disturbance_years_csv.jl")
include("disturbance_data_gen/02_a_create_cover_manta_csv.jl")
include("disturbance_data_gen/02_b_create_cover_transect_csv.jl")

include("disturbance_data_gen/03_fill_missing_transect_disturbances.jl")
include("disturbance_data_gen/04_create_survival_rates_csv.jl")
include("disturbance_data_gen/05_build_dist_mortality_datacube.jl")