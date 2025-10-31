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

# Create template df to fill with disturbance start and end years
# Only disturbances marked as storms or cots or multiple (in the last case, only if
# the description mentions "storm" or "cots") are considered in this dataset

# Manta tow template csv
_disturbance_years_df_manta = create_template_disturbance_years(reef_name_and_ids, ["m", "s", "c"])

# Remove rows that are multiple but don't mention storm or cots in their description
disturbances_description_manta = lowercase.(_disturbance_years_df_manta.description)
has_storm_description_manta = occursin.(r"yc", disturbances_description_manta) .||
                              occursin.(r"storm", disturbances_description_manta)
has_cots_description_manta = occursin.(r"cots", disturbances_description_manta) .||
                             occursin.(r"crown", disturbances_description_manta)
is_multiple_disturbance_manta = _disturbance_years_df_manta.type .== "m"

target_remove_rows_manta = is_multiple_disturbance_manta .&& (.!(has_storm_description_manta .|| has_cots_description_manta))
if sum(target_remove_rows_manta) > 0
    deleteat!(_disturbance_years_df_manta, findfirst(target_remove_rows_manta))
end

# Save disturbances_df as template csv
disturbance_years_manta_df = DataFrames.select(sort(_disturbance_years_df_manta, :reef_name), Not(:description))

# The reason why this file has a "template" in the name is: after I manually fill it I rename
# it to remove the word "template" so that, if I run this script again accidentally,
# I don't accidentally lose all the data I've manually collected (:
CSV.write("data/template_disturbance_years_manta.csv", disturbance_years_manta_df)

# Transect template csv
_disturbance_years_transect_df = create_template_disturbance_years(reef_name_and_ids, ["m", "s", "c"]; sample_type="PPOINT")

# Remove rows that are multiple but don't mention storm or cots in their description
disturbances_description_transect = lowercase.(_disturbance_years_transect_df.description)
has_storm_description_transect = occursin.(r"yc", disturbances_description_transect) .||
                                 occursin.(r"storm", disturbances_description_transect)
has_cots_description_transect = occursin.(r"cots", disturbances_description_transect) .||
                                occursin.(r"crown", disturbances_description_transect)
is_multiple_disturbance_transect = _disturbance_years_transect_df.type .== "m"

target_remove_rows_transect = is_multiple_disturbance_transect .&& (.!(has_storm_description_transect .|| has_cots_description_transect))
if sum(target_remove_rows_transect) > 0
    deleteat!(_disturbance_years_transect_df, findfirst(target_remove_rows_transect))
end

# Save disturbances_df as template csv
disturbance_years_transect_df = DataFrames.select(sort(_disturbance_years_transect_df, :reef_name), Not(:description))
CSV.write("data/template_disturbance_years_transect.csv", disturbance_years_transect_df)