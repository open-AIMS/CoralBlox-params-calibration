# `rm` = ReefMonitoring. Some target_loc_ids don't match any ReefMonitoring reef name -
# those rows get an empty reef_names entry and are dropped later.
rm_reef_names, rm_reef_ids = rm_reef_spec(canonical_gpkg)
tmp_target_loc_names = [rm_reef_names[rm_reef_ids.==id] for id in target_loc_ids]
target_loc_names = [isempty(name) ? "" : name[1] for name in tmp_target_loc_names]
reef_name_and_ids = DataFrame((reef_names=target_loc_names, reef_ids=target_loc_ids))

# Only storms, COTS, or "multiple" disturbances whose description mentions storm/cots.
_disturbance_years_df_manta = create_template_disturbance_years(reef_name_and_ids, ["m", "s", "c"])

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

disturbance_years_manta_df = DataFrames.select(sort(_disturbance_years_df_manta, :reef_name), Not(:description))

# "template" in the filename: after manually filling year_before/year_after, the curated copy
# is renamed to drop "template" - so rerunning this script never clobbers curated data.
CSV.write(joinpath(DATA_DIR, "template_disturbance_years_manta.csv"), disturbance_years_manta_df)

_disturbance_years_transect_df = create_template_disturbance_years(reef_name_and_ids, ["m", "s", "c"]; sample_type="PPOINT")
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

disturbance_years_transect_df = DataFrames.select(sort(_disturbance_years_transect_df, :reef_name), Not(:description))
CSV.write(joinpath(DATA_DIR, "template_disturbance_years_transect.csv"), disturbance_years_transect_df)