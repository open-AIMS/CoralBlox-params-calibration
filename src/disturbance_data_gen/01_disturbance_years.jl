using ReefMonitoring
using DataFrames, CSV

include("base.jl")

# Names and ids for all reefs on Reef Monitoring API. `rm` here stands for ReefMonitoring
rm_reef_names, rm_reef_ids = rm_reef_spec(canonical_gpkg)

# Get reef names for `target_loc_id` reefs.
# Some `target_loc_id`s don't match any of the ReefMonitoring ids
tmp_target_loc_names = [rm_reef_names[rm_reef_ids.==id] for id in target_loc_ids]
target_loc_names = [isempty(name) ? "" : name[1] for name in tmp_target_loc_names]

# Dataframe to match each `target_loc_id` with a `rm_reef_name`
reef_name_and_ids = DataFrame((reef_names=target_loc_names, reef_ids=target_loc_ids))

# Get the ids of the target reefs that didn't match any rm_reef_name:
# no_match_ids = reef_name_and_ids[reef_name_and_ids.reef_names.=="", :].reef_ids
# ["14116101104", "21278100104"]

# TODO I used the code below twice to generate manta tow and transect files. Needs improvement
# Create template df to fill with disturbance start and end years
# Only manta disturbances marked as storms or cots or multiple (in the last case, only if
# the description mentions "storm" or "cots") are considered in this dataset
disturbance_years_df = create_template_disturbance_years(reef_name_and_ids, ["m", "s", "c"])

# Remove rows that are multiple but don't mention storm or cots in their description
disturbances_description = lowercase.(disturbance_years_df.description)
has_storm_description = occursin.(r"yc", disturbances_description) .||
                        occursin.(r"storm", disturbances_description)
has_cots_description = occursin.(r"cots", disturbances_description) .||
                       occursin.(r"crown", disturbances_description)
is_multiple_disturbance = disturbance_years_df.type .== "m"

target_remove_rows = is_multiple_disturbance .&& (.!(has_storm_description .|| has_cots_description))
if sum(target_remove_rows) > 0
    deleteat!(disturbance_years_df, findfirst(target_remove_rows))
end

# Save disturbances_df as template csv
manta_disturbance_years_df = DataFrames.select(sort(disturbance_years_df, :reef_name), Not(:description))

# The reason why this file has a "template" in the name is: after I manually fill it I rename
# it to remove the word "template" so that, if I run this script again accidentally,
# I don't accidentally lose all the data I've manually collected (:
CSV.write("template_disturbance_years_manta.csv", manta_disturbance_years_df)
