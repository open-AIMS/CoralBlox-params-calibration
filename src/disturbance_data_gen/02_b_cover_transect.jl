using YAXArrays, NetCDF
using CSV, DataFrames

# After running stript 01, the "template_disturbance_years_transect.csv"'s columns
# `year_before` and `year_after` of the , were manually filled by eye picking at the
# Reef Monitoring Dashboard (https://apps.aims.gov.au/reef-monitoring/reefs)
# and put into the file "disturbance_years_transect.csv"

# Add `survival_rate`s to `disturbance_years_transect.csv`
disturbance_years_transect_path = "data/disturbance_years_transect.csv"
disturbance_years_transect = CSV.read(disturbance_years_transect_path, DataFrame; stringtype=String, comment="#")

# Value in COMPOSITION_PATH were extracted from ReefMonitoring API and are used to fill the
# benthic composition to disturbance_years_transect
coral_composition = open_dataset(COMPOSITION_PATH).mean
composition_reef_ids = coral_composition.location.val.data

functional_groups = coral_composition.taxa.val.data
for (group) in functional_groups
    disturbance_years_transect[!, "cover_before_$(group)"] .= 0.0
    disturbance_years_transect[!, "cover_after_$(group)"] .= 0.0
end

transect_reefs = unique(disturbance_years_transect[:, [:reef_name, :reef_id]])
transect_reef_names = transect_reefs.reef_name
transect_reef_ids = string.(transect_reefs.reef_id)

cache_mask = falses(nrow(disturbance_years_transect))
for (reef_idx, reef_name) in enumerate(transect_reef_names)
    reef_id = transect_reef_ids[reef_idx]
    if !(reef_id in composition_reef_ids)
        @info("Reef ID not found in coral_composition: $reef_id")
        continue
    end

    target_cover = coral_composition[location=At(reef_id), taxa=At(functional_groups)]
    #transect_cover = ReefMonitoring.get_photo_transect(reef_name)

    # Find cyclones for current reef
    cache_mask .= disturbance_years_transect.reef_name .== reef_name
    reef_disturbances = disturbance_years_transect[cache_mask, :]

    # If there are two disturbances in the same year, the for loop below will not work
    @assert length(reef_disturbances.year) == length(unique(reef_disturbances.year))

    for disturbance in eachrow(reef_disturbances)
        year_before, year_after = disturbance.year_before, disturbance.year_after

        # Mask to match by reef name and disturbance year
        cache_mask .= disturbance_years_transect.reef_name .== reef_name
        cache_mask .= cache_mask .&& (disturbance_years_transect.year .== disturbance.year)

        try
            cover_before = read(target_cover[timesteps=At(year_before)])
            cover_after = read(target_cover[timesteps=At(year_after)])
            disturbance_years_transect[cache_mask, cover_before_cols] .= cover_before'
            disturbance_years_transect[cache_mask, cover_after_cols] .= cover_after'
        catch
            @info("Reef name: $reef_name \nCyclone: $cyclone")
            continue
        end
    end
end

# Assert there are no NaN values
@assert all(Matrix(.!isnan.(disturbance_years_transect[!, cover_before_cols])))
@assert all(Matrix(.!isnan.(disturbance_years_transect[!, cover_after_cols])))

# Assert there are no survival rate greater than 1
@assert all((0.0 .<= Matrix(disturbance_years_transect[!, cover_before_cols]) .<= 1.0))
@assert all((0.0 .<= Matrix(disturbance_years_transect[!, cover_after_cols]) .<= 1.0))

# Manually split cover for these two reefs in two steps
split_taxa_cover!(disturbance_years_transect, "Taylor Reef", 2016, 2018)
split_taxa_cover!(disturbance_years_transect, "Hoskyn Island", 2008, 2010)

# Delete rows with all zero values
zero_before_mask = sum(Matrix(disturbance_years_transect[!, cover_before_cols]), dims=2) .== 0.0
zero_after_mask = sum(Matrix(disturbance_years_transect[!, cover_after_cols]), dims=2) .== 0.0
deleteat!(disturbance_years_transect, getindex.(findall(zero_before_mask .&& zero_after_mask), 1))

CSV.write("historic_cyclone_trajectories/disturbance_cover_transect.csv", disturbance_years_transect)
