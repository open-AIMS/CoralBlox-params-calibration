using ADRIA
using CSV, YAXArrays, NetCDF, DataFrames
using Dates, TimeZones

survival_rates_path = joinpath(DATA_DIR, "disturbance_survival_rates.csv")
survival_rates_df = CSV.read(survival_rates_path, DataFrame; stringtype=String, comment="#")

# Ensure that there are no duplicate pairs of reef/year
@assert length(unique(eachrow(survival_rates_df[:, [:reef_name, :dist_year]]))) == nrow(survival_rates_df)

START_YEAR, END_YEAR = 2008, 2022

dom = ADRIA.load_domain(RMEDomain, RME_DOMAIN_PATH, "45", timeframe=(START_YEAR, END_YEAR))
new_cyclone_mortality_scens = deepcopy(dom.cyclone_mortality_scens)

# Convert RME_UNIQUE_ID to RME_GBRMPA_ID because that's what is used in the RME domain
target_gbrmpa_ids = canonical_gpkg[[findfirst(x -> x == id, canonical_gpkg.RME_UNIQUE_ID) for id in target_loc_ids], :RME_GBRMPA_ID]

# Set disturbance mortality rate to zero for all calibration locations
new_cyclone_mortality_scens[locs=At(target_gbrmpa_ids)] .= 0.0

# disturbance_reef_ids are not the same as target_loc_ids. The former is a subset of the latter.
disturbance_reef_ids = unique(survival_rates_df.reef_id)

survival_cols = "survival_rate_" .* functional_groups

for survival_col in survival_cols
    survival_rates_df[:, survival_col] .= replace!(survival_rates_df[:, survival_col], NaN => 0.0)
end

for reef_id in disturbance_reef_ids
    survival_df = survival_rates_df[survival_rates_df.reef_id.==reef_id, vcat(:dist_year, Symbol.(survival_cols))]

    reef_gbrmpa_id = canonical_gpkg[canonical_gpkg.RME_UNIQUE_ID.==string.(reef_id), :RME_GBRMPA_ID][1]
    for survival_rate in eachrow(survival_df)
        year = max(Int64(survival_rate.dist_year), 2008)

        #! Maximum allowed mortality rate is 0.95
        mortality_rates = clamp.(1 .- collect(survival_rate[Symbol.(survival_cols)]), 0.0, 0.95)
        new_cyclone_mortality_scens[locs=At(reef_gbrmpa_id), timesteps=At(year)] .= mortality_rates
    end
end

new_axes = (
    Dim{:timesteps}(2008:2022),
    Dim{:locs}(collect(new_cyclone_mortality_scens.locs)),
    Dim{:species}(String.([
        :tabular_Acropora,
        :corymbose_Acropora,
        :corymbose_non_Acropora,
        :small_massives,
        :large_massives
    ])),
    Dim{:scenarios}(1:2)
)

properties = Dict{String,String}(
    "timesteps" => "Vector{Int64} 2008:2022",
    "locs" => "Vector{String} ['10-330', …, '23-049']",
    "species" => "Vector{String} ['tabular_Acropora', 'corymbose_Acropora', 'corymbose_non_Acropora', 'small_massives', 'large_massives']",
    "scenarios" => "Vector{Int64} [1,2] (both scenarios are identical)",
    "unit" => "Proportion of coral cover lost (0-1)",
    "created_at" => string(now(tz"UTC")),
)

disturbance_mortality_scens = YAXArray(
    new_axes,
    repeat(new_cyclone_mortality_scens.data, 1, 1, 1, 2),
    properties
)

savedataset(
    Dataset(; disturbance_mortality_scens);
    path=joinpath(DATA_DIR, "historical_disturbance_mortality_rates", "historical_disturbance_mortality_rates.nc"),
    driver=:netcdf,
    overwrite=true
)
