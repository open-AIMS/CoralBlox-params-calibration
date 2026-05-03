using Revise, Infiltrator
using Serialization

using TOML
using Dates
using CSV, DataFrames, YAXArrays
using StatsBase, Statistics

using ADRIA
using ADRIA: GDF, AG, DimensionalData

src_path = dirname(@__DIR__)
datasets_path = joinpath(dirname(src_path), "datasets")

global CONFIG = TOML.parsefile(joinpath(src_path, "calib_config.toml"))

# Configuration sections
global DOMAIN_CONFIG = CONFIG["Domains"]
global GEOSPATIAL_CONFIG = CONFIG["Geospatial"]
global INITIALISATION_CONFIG = CONFIG["Initialisation"]
global OUTPUT_CONFIG = CONFIG["Outputs"]

# ADRIA Domain paths
global RME_DOMAIN_PATH = DOMAIN_CONFIG["rme_domain"]
global HISTORICAL_DHW_PATH = joinpath(
    src_path,
    "historical_dhw_data_gen",
    "data",
    "historical_dhw.nc"
)
global HISTORICAL_CYCLONE_MORTALITY_PATH = joinpath(
    src_path,
    "historical_disturbance_mortality_data_gen/",
    "data/",
    "historical_disturbance_mortality_rates/",
    "historical_disturbance_mortality_rates.nc"
)

# Geospatial filepaths
global CANONICAL_PATH = GEOSPATIAL_CONFIG["canonical_path"]
global LTMP_SHP_PATH = get(
    GEOSPATIAL_CONFIG,
    "ltmp_shp",
    joinpath(datasets_path, "spatial_data/gbr_3Zone 2.shp")
)
global LOC_CLASS_PATH = get(
    GEOSPATIAL_CONFIG,
    "classification_path",
    joinpath(datasets_path, "spatial_data/location_classification_MPA.csv")
)

global GBRMPA_MAINLAND_PATH = joinpath(datasets_path, "spatial_data/Great_Barrier_Reef_Features.geojson")

# Calibration Target / Observational Data
global TARGET_CONFIG = get(CONFIG, "Observations", Dict())
global LOC_CLASS_TARGET_PATH = get(
    TARGET_CONFIG,
    "manta_tow_path",
    joinpath(datasets_path, "ltmp_data/manta_tow_mean_std.nc")
)
global LTMP_REEF_DATA_PATH = get(
    TARGET_CONFIG,
    "ltmp_reef_data",
    joinpath(datasets_path, "ltmp_data/manta_tow_data_reef_lvl.gpkg")
)
global COMPOSITION_PATH = get(
    TARGET_CONFIG,
    "composition_netcdf",
    joinpath(datasets_path, "ltmp_data/coral_composition.nc")
)
global LTMP_MODELLED_OBS_PATH = get(
    TARGET_CONFIG,
    "ltmp_modelled_obs",
    joinpath(datasets_path, "ltmp_data/modelled_brms.beta.ry.disp.csv")
)

# Initialisation filepaths
global INIT_COVER_PATH = get(
    INITIALISATION_CONFIG,
    "init_cover_filepath",
    joinpath(datasets_path, "spatial_data/init_cover.dat")
)
global INIT_GUESS_PATH = get(INITIALISATION_CONFIG, "init_guess_filepath", "")

# Output filepaths
global OUT_DIR = OUTPUT_CONFIG["out_dir"]
global RESULT_FN = get(OUTPUT_CONFIG, "result_filename", "results.dat")

const START_YEAR = 2008
const END_YEAR = 2022

if !@isdefined(OPTIONS)
    global reload_domain = false

    # define OPTIONS to prevent reinitialise on every include
    global OPTIONS = true
end

# julia> include("path to gist.jl")
# julia> create_compare_plot("plot file name.png")
# julia> create_error_statistics("error statistics filename.csv")
# julia> plot_residuals("residual plot.png")

# Convenience functions to retrieve parameters from an ADRIA model specification
function extract_param_group_idx(model_spec::DataFrame, needle::String)::Vector{Int64}
    needle_pos = contains.(string.(model_spec.fieldname), needle)
    return findall(needle_pos)
end
function extract_param_group_idx(
    model::ADRIA.Model, component, needle::String
)::Vector{Int64}
    comp_params = ADRIA.component_params(model, component)
    return extract_param_group_idx(comp_params, needle)
end

function extract_param_group(model::ADRIA.Model, component, needle::String)::DataFrame
    group_pos = extract_param_group_idx(model, component, needle)
    comp_params = ADRIA.component_params(model, component)
    return comp_params[group_pos, :]
end
function extract_param_group(model_spec::DataFrame, needle::String)::DataFrame
    group_pos = extract_param_group_idx(model_spec, needle)
    return model_spec[group_pos, :]
end

"""
To minimize the number of index vectors being passed around and copies being created, define
an immutable struct to provide an interface to move between indexing between different
arrays the refer to the same locations.
"""
struct LocationDataStore
    # Data fields
    domain_gpkg::DataFrame
    ltmp_unique_ids::Vector{String}
    ltmp_coral_cover::Array{Union{Missing,Float64},2}
    coral_composition::Array{Union{Missing,Float64},3}
    # Index Fields
    ltmp_cover_to_domain::Vector{Int64}
    composition_to_domain::Vector{Int64}
end

"""
    ltmp_cover_idx_to_domain(loc_data_store::LocationDataStore, ltmp_cover_idx::Int64)

Given an index of a location in the ltmp reef coral cover dataframe, return the index of the
same location in the domain geopackage.
"""
function ltmp_cover_idx_to_domain(loc_data_store::LocationDataStore, ltmp_cover_idx::Int64)
    return loc_data_store.ltmp_cover_to_domain[ltmp_cover_idx]
end

"""
    composition_idx_to_domain(loc_data_store::LocationDataStore, composition_idx::Int64)

Given an index of a location in the coral composition yaxarray, return the index of the
same location in the domain geopackage.
"""
function composition_idx_to_domain(
    loc_data_store::LocationDataStore, composition_idx::Int64
)
    return loc_data_store.composition_to_domain[composition_idx]
end

function get_ltmp_loc_unique_id(loc_data_store::LocationDataStore, ltmp_idx::Int64)::String
    return loc_data_store.ltmp_unique_ids[ltmp_idx]
end

function get_composition_loc_unique_id(
    loc_data_store::LocationDataStore, composition_idx::Int64
)::String
    return loc_data_store.coral_composition.location[composition_idx]
end

include("./perf_metrics.jl")
