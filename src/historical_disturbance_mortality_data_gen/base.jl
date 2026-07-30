using ADRIA
using ReefMonitoring
import GeoDataFrames as GDF
using CSV, DataFrames
using YAXArrays, NetCDF
using TOML

isdefined(Main, :DATA_DIR) || (global DATA_DIR = joinpath(@__DIR__, "data"))

isdefined(Main, :CONFIG) || (global CONFIG = TOML.parsefile(joinpath(@__DIR__, "..", "..", "config.toml")))

isdefined(Main, :OUTPUT_CONFIG) || (global OUTPUT_CONFIG = CONFIG["calibration"]["outputs"])
isdefined(Main, :OUT_DIR) || (global OUT_DIR = OUTPUT_CONFIG["out_dir"])

isdefined(Main, :DOMAIN_CONFIG) || global DOMAIN_CONFIG = CONFIG["calibration"]["domains"]
isdefined(Main, :RME_DOMAIN_PATH) || global RME_DOMAIN_PATH = DOMAIN_CONFIG["rme_domain"]

isdefined(Main, :GEOSPATIAL_CONFIG) || (global GEOSPATIAL_CONFIG = CONFIG["calibration"]["geospatial"])
isdefined(Main, :CANONICAL_PATH) || (global CANONICAL_PATH = GEOSPATIAL_CONFIG["canonical_path"])

# Kept in sync manually with CalibrationConfig.composition_path/ltmp_reef_data_path in
# src/common/calib_setup.jl - this pipeline predates that config API.
isdefined(Main, :TARGET_CONFIG) || (global TARGET_CONFIG = get(CONFIG["calibration"], "observations", Dict()))
isdefined(Main, :COMPOSITION_PATH) || (global COMPOSITION_PATH = get(
    TARGET_CONFIG, "composition_netcdf",
    joinpath(@__DIR__, "..", "..", "datasets", "ltmp_data", "coral_composition.nc")
))
isdefined(Main, :LTMP_REEF_DATA_PATH) || (global LTMP_REEF_DATA_PATH = get(
    TARGET_CONFIG, "ltmp_reef_data",
    joinpath(@__DIR__, "..", "..", "datasets", "ltmp_data", "manta_tow_data_reef_lvl.gpkg")
))

"""
Ids of every calibration+validation location - every reef mortality needs to be computed for.
"""
function _target_loc_ids(calib_split_path)::Vector{String}
    calib_split = CSV.read(calib_split_path, DataFrame)
    return string.(sort(calib_split.UNIQUE_IDS))
end

"""
Ids of CALIBRATION-only locations. Used for growth-rate estimation, so validation reefs are
never used to fit the correction their own mortality is later evaluated against.
"""
function _calibration_loc_ids(calib_split_path)::Vector{String}
    calib_split = CSV.read(calib_split_path, DataFrame)
    calibration_only = calib_split[calib_split.USAGE.=="calibration", :]
    return string.(sort(calibration_only.UNIQUE_IDS))
end

# Cached locally rather than under config.toml's out_dir, to avoid a circular dependency on
# historical_disturbance_mortality_rates.nc (the file this pipeline builds). See
# base/generate_calibration_split.jl.
local_calib_split_path = joinpath(DATA_DIR, "calibration_split.csv")
if isfile(local_calib_split_path)
    @info "Found local_calib_split file, skipping its generation."
else
    @info "Generating local_calib_split file."
    include(joinpath(@__DIR__, "base", "generate_calibration_split.jl"))
    @info "local_calib_split file successfully generated."
end

target_loc_ids = _target_loc_ids(local_calib_split_path)
calibration_loc_ids = _calibration_loc_ids(local_calib_split_path)

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

# Caches ReefMonitoring.get_manta_tow's raw output, shared by 03_a and
# generate_clean_growth_intervals.jl below. Delete to force a fresh pull.
manta_tow_cache_path = joinpath(DATA_DIR, "manta_tow_raw_cache.csv")
if isfile(manta_tow_cache_path)
    @info "Found manta_tow_raw_cache file, skipping its generation."
else
    @info "Generating manta_tow_raw_cache file."
    include(joinpath(@__DIR__, "base", "generate_manta_tow_cache.jl"))
    @info "manta_tow_raw_cache file successfully generated."
end

# Guards on background_growth_rates.csv, not clean_growth_intervals.csv's own existence -
# clean_growth_intervals.csv is only ever a means to the former. See
# base/generate_clean_growth_intervals.jl.
background_growth_rates_path = joinpath(DATA_DIR, "background_growth_rates.csv")
if isfile(background_growth_rates_path)
    @info "Found background_growth_rates file, its skipping generation."
else
    @info "Generating background_growth_rates file."
    include(joinpath(@__DIR__, "base", "generate_clean_growth_intervals.jl"))
    @info "background_growth_rates file successfully generated."
end

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
Splits a manta cover interval that spans two undetected disturbances into two, inserting
`intermediate_year` as a hook to attach the missing disturbance's data to.

```julia
split_cover!(disturbance_manta_years, "Rebe Reef", 2005, 2011, 2009)
```
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

        n_target_disturbances = nrow(target_disturbances)
        storm_names = replace(target_disturbances.storm_name, nothing => "unknown")
        description = replace(target_disturbances.description, nothing => "")
        storm_years = parse.(Float64, target_disturbances.ddate)

        new_entries = hcat(
            storm_names,
            target_disturbances.disturbance,
            storm_years,
            description,
            repeat([reef_name reef_id 0 0],
                n_target_disturbances)
        )

        push!.(Ref(template_disturbance_years_df), eachrow(new_entries))
    end

    return template_disturbance_years_df
end
