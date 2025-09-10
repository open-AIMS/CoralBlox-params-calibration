using DataFrames, CSV
using NetCDF
import YAXArrays as YA
using ReefMonitoring

# include("common.jl")

# After running stript 01, the "template_disturbance_years_manta.csv"'s columns
# `year_before` and `year_after` of the , were manually filled by eye picking at the
# Reef Monitoring Dashboard (https://apps.aims.gov.au/reef-monitoring/reefs)
# and put into the file "disturbance_years_manta.csv"

# Add `survival_rate`s to `disturbance_years_manta.csv`
distribution_manta_years_path = "data/disturbance_years_manta.csv"
disturbance_manta_years = CSV.read(distribution_manta_years_path, DataFrame; stringtype=String, comment="#")
disturbance_manta_years.cover_before = zeros(Float64, nrow(disturbance_manta_years))
disturbance_manta_years.cover_after = zeros(Float64, nrow(disturbance_manta_years))

cache_mask = falses(nrow(disturbance_manta_years))
for reef_name in unique(disturbance_manta_years.reef_name)
    # Use Reef Monitoring API to get manta tow cover data
    manta_cover = ReefMonitoring.get_manta_tow(reef_name)

    # Find cyclones for current reef
    cache_mask .= disturbance_manta_years.reef_name .== reef_name
    reef_cyclones = disturbance_manta_years[cache_mask, :]

    for cyclone in eachrow(reef_cyclones)
        year_before, year_after = cyclone.year_before, cyclone.year_after

        cache_mask .= disturbance_manta_years.reef_name .== reef_name
        cache_mask .= cache_mask .&& (disturbance_manta_years.cyclone_name .== cyclone.cyclone_name)

        try
            disturbance_manta_years[cache_mask, :cover_before] .= manta_cover[manta_cover.report_year.==year_before, :mean][1]
            disturbance_manta_years[cache_mask, :cover_after] .= manta_cover[manta_cover.report_year.==year_after, :mean][1]
        catch
            @info("Reef name: $reef_name \nCyclone: $cyclone")
            continue
        end
    end
end

# Manually add data for Rebe Reef cyclone Hamish
# Add a "fake" data point for reefs with large time gaps in manta data to account for additional disturbances.
split_cover!(disturbance_manta_years, "Rebe Reef", 2005, 2011, 2009)
split_cover!(disturbance_manta_years, "Martin Reef", 2013, 2017, 2015)

CSV.write("data/disturbance_cover_manta.csv", disturbance_manta_years)
