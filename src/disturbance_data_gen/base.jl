using ADRIA
using ReefMonitoring
import GeoDataFrames as GDF
using CSV, DataFrames
using YAXArrays, NetCDF
using TOML

isdefined(Main, :CONFIG) || (global CONFIG = TOML.parsefile("calib_config.toml"))

isdefined(Main, :OUTPUT_CONFIG) || (global OUTPUT_CONFIG = CONFIG["Outputs"])
isdefined(Main, :OUT_DIR) || (global OUT_DIR = OUTPUT_CONFIG["out_dir"])

isdefined(Main, :DOMAIN_CONFIG) || global DOMAIN_CONFIG = CONFIG["Domains"]
isdefined(Main, :RME_DOMAIN_PATH) || global RME_DOMAIN_PATH = DOMAIN_CONFIG["rme_domain"]

isdefined(Main, :GEOSPATIAL_CONFIG) || (global GEOSPATIAL_CONFIG = CONFIG["Geospatial"])
isdefined(Main, :CANONICAL_PATH) || (global CANONICAL_PATH = GEOSPATIAL_CONFIG["canonical_path"])

"""
Ids of the locations we use in the calibration
"""
function _target_loc_ids(calib_split_path)::Vector{String}
    calib_split = CSV.read(calib_split_path, DataFrame)
    return string.(sort(calib_split.UNIQUE_IDS))
end
target_loc_ids = _target_loc_ids(joinpath(OUT_DIR, "calibration_split.csv"))

if !@isdefined(canonical_gpkg)
    @info "Loading Canonical gpkg"
    canonical_gpkg = GDF.read(CANONICAL_PATH)
end

# Could get these from ADRIA.functional_groups but don't want to load ADRIA just for this
functional_groups = [
    "Tabular_Acropora",
    "Corymbose_Acropora",
    "Corymbose_non_Acropora",
    "Small_Massives",
    "Large_Massives"
]

cover_before_col_names = "cover_before_" .* functional_groups
cover_after_col_names = "cover_after_" .* functional_groups

"""
    rm_reef_spec(canonical_gpkg)::Tuple{Vector{String},Vector{String}}

Names and ids for all reefs on Reef Monitoring API
"""
function rm_reef_spec(canonical_gpkg)::Tuple{Vector{String},Vector{String}}
    spec = ReefMonitoring.get_reef_info()
    reef_ids::Vector{String} = ReefMonitoring.get_unique_id.(
        spec.latitude,
        spec.longitude,
        [canonical_gpkg]
    )
    return spec.aims_reef_name, reef_ids
end

"""
For some reefs there are disturbances reported in the transect data that are not reported
in the manta data. In some of these cases, there is a big time window in the manta without
data (e.g. in Rebe Reef there is one data point in 2005 and another in 2011) with a big
drop in cover caused by two disturbances. For these, I add an extra datapoint to have a
"hook" where I can add the data for this extra disturbance.

# Example
```julia
split_cover!(disturbance_manta_years, "Rebe Reef", 2005, 2011, 2009)
````
"""
function split_cover!(
    disturbance_cover,
    reef_name,
    year_before,
    year_after,
    intermediate_year,
)::Nothing
    target_mask = (disturbance_cover.reef_name .== reef_name) .&&
                  (disturbance_cover.year_before .== year_before) .&&
                  (disturbance_cover.year_after .== year_after)
    tmp_df = disturbance_cover[target_mask, :]

    split_year_first, split_year_last = extrema(tmp_df.year)

    disturbance_first = tmp_df[tmp_df.year.==split_year_first, :]
    disturbance_last = tmp_df[tmp_df.year.==split_year_last, :]

    cover_first_before = disturbance_first.cover_before[1]
    cover_last_after = disturbance_last.cover_after[1]
    cover_interemediate = (cover_first_before + cover_last_after) / 2

    mask_first = target_mask .&& disturbance_cover.year .== disturbance_first.year
    mask_last = target_mask .&& disturbance_cover.year .== disturbance_last.year
    disturbance_cover[mask_first, :cover_after] .= cover_interemediate
    disturbance_cover[mask_last, :cover_before] .= cover_interemediate

    disturbance_cover[mask_first, :year_after] .= intermediate_year
    disturbance_cover[mask_last, :year_before] .= intermediate_year

    return nothing
end

function split_taxa_cover!(disturbance_cover, reef_name, year_before, year_after)::Nothing
    target_mask = (disturbance_cover.reef_name .== reef_name) .&&
                  (disturbance_cover.year_before .== year_before) .&&
                  (disturbance_cover.year_after .== year_after)
    tmp_df = disturbance_cover[target_mask, :]

    split_year_low, split_year_up = extrema(tmp_df.year)

    disturbance_low = tmp_df[tmp_df.year.==split_year_low, :]
    disturbance_up = tmp_df[tmp_df.year.==split_year_up, :]

    cover_before_low = Matrix(disturbance_low[:, Symbol.(cover_before_col_names)])
    cover_after_up = Matrix(disturbance_up[:, Symbol.(cover_after_col_names)])
    cover_interemediate = (cover_before_low .+ cover_after_up) ./ 2

    mask_low = target_mask .&& disturbance_cover.year .== disturbance_low.year
    mask_up = target_mask .&& disturbance_cover.year .== disturbance_up.year
    disturbance_cover[mask_low, cover_after_col_names] .= cover_interemediate
    disturbance_cover[mask_up, cover_before_col_names] .= cover_interemediate

    return nothing
end

"""
    create_template_disturbance_years(
        reef_name_and_ids::DataFrame, disturbance_type::Vector{String};
        sample_type::String="MANTA"
    )::DataFrame

# Arguments
- `reef_name_and_ids::DataFrame` :
- `disturbance_type` : Valid options: `"s"` for storms, `"m"` for multiple, etc
- `sample_type::String="MANTA"` : Valid options: `"MANTA"`, `"PPOINT"`
"""
function create_template_disturbance_years(
    reef_name_and_ids::DataFrame, disturbance_type::Vector{String};
    sample_type::String="MANTA"
)::DataFrame
    template_disturbance_years_df = DataFrame(
        cyclone_name=String[],
        type=String[],
        year=Float64[],
        description=String[],
        reef_name=String[],
        reef_id=String[],
        year_before=Int64[],        # empty, to be filled, year before the disturbance
        year_after=Int64[]          # empty, to be filled, year after the disturbance
    )

    start_year = 2008
    end_year = 2022

    # Use the ReefMonitoring api to get disturbance data for each reef considered (excluding
    # reefs for which there is no reef_name associated with that reef's reef_id)
    # For each disturbance, adds one entry to template_disturbance_years_df
    # Each entry will become a row on the template CSV
    for reef_row in eachrow(reef_name_and_ids)
        reef_name, reef_id = reef_row
        isempty(reef_name) ? continue : 0

        disturbances::Union{DataFrame,Nothing} = ReefMonitoring.get_disturbances(reef_name)
        isempty(disturbances) ? continue : 0

        is_target_sample_type = disturbances.sample_type .== sample_type
        is_target_disturbance_type = disturbances.disturbance .∈ [disturbance_type]
        is_within_time_interval = start_year .<= parse.(Float64, disturbances.ddate) .<= end_year

        target_disturbances =
            disturbances[
                is_target_sample_type.&&is_target_disturbance_type.&&is_within_time_interval,
                :
            ]

        # Build new `template_disturbance_years_df` entries
        n_target_disturbances = nrow(target_disturbances)
        storm_names = replace(target_disturbances.storm_name, nothing => "unknown")
        description = replace(target_disturbances.description, nothing => "")
        storm_years = parse.(Float64, target_disturbances.ddate)

        # Collate results to build new entries
        new_entries = hcat(
            storm_names,
            target_disturbances.disturbance,
            storm_years,
            description,
            repeat([reef_name reef_id 0 0],
                n_target_disturbances)
        )

        # Add new entries to `template_disturbance_years_df`
        push!.(Ref(template_disturbance_years_df), eachrow(new_entries))
    end

    return template_disturbance_years_df
end