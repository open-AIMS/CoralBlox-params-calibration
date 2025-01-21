using Revise, Infiltrator
using Serialization

using TOML
using Dates
using CSV, DataFrames, YAXArrays
using StatsBase, Statistics
using BlackBoxOptim

using CairoMakie
using ADRIA
using ADRIA: GDF, AG, DimensionalData

function ADRIA.bin_edges(; unit=:m)
    return Matrix(
        [
            2.5 7.5 12.5 25.0 50.0 80.0 120.0 160.0;
            2.5 7.5 12.5 20.0 30.0 60.0 100.0 150.0;
            2.5 7.5 12.5 20.0 30.0 40.0 50.0 60.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0
        ]
    ) .* ADRIA.linear_scale(:cm, unit)
end


global CONFIG = TOML.parsefile("calib_config.toml")

# Configuration sections
global DOMAIN_CONFIG         = CONFIG["Domains"]
global GEOSPATIAL_CONFIG     = CONFIG["Geospatial"]
global TARGET_CONFIG         = CONFIG["Observations"]
global INITIALISATION_CONFIG = CONFIG["Initialisation"]
global OUTPUT_CONFIG         = CONFIG["Outputs"]

# ADRIA Domain paths
global REEFMOD_DOMAIN_PATH = DOMAIN_CONFIG["reefmod_domain"]
global RME_DOMAIN_PATH     = DOMAIN_CONFIG["rme_domain"]
global HISTORIC_DHW_PATH   = joinpath(RME_DOMAIN_PATH, "data_files", "dhw", "GBR_past_DHW_CRW_5km_1985_2022_Dec_2022.csv")

# Geospatial filepaths
global CANONICAL_PATH        = GEOSPATIAL_CONFIG["canonical_path"]
global LTMP_SHP_PATH         = GEOSPATIAL_CONFIG["ltmp_shp"]
global LOC_CLASS_PATH        = GEOSPATIAL_CONFIG["classification_path"] # location classes
global BIOREGION_GROUPS_PATH = GEOSPATIAL_CONFIG["bioregion_group_gpkg"]

# Calibration Target / Observational Data
global LOC_CLASS_TARGET_PATH  = TARGET_CONFIG["manta_tow_path"]  # target data for location classes
global LTMP_REEF_DATA_PATH    = TARGET_CONFIG["ltmp_reef_data"]  # target data for ltmp locs
global COMPOSITION_PATH       = TARGET_CONFIG["composition_netcdf"]
global LTMP_MODELLED_OBS_PATH = TARGET_CONFIG["ltmp_modelled_obs"]

# Initialisation filepaths
global INIT_COVER_PATH    = INITIALISATION_CONFIG["init_cover_filepath"]
global INIT_GUESS_PATH    = INITIALISATION_CONFIG["init_guess_filepath"]  # optional
global ECORRAP_PARAM_PATH = INITIALISATION_CONFIG["ecorrap_param_filepath"]

# Output filepaths
global OUT_DIR   = OUTPUT_CONFIG["out_dir"]
global RESULT_FN = OUTPUT_CONFIG["result_filename"]

const START_YEAR = 2008
const END_YEAR = 2022
global n_scenarios = 16

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
function extract_param_group_idx(model::ADRIA.Model, component, needle::String)::Vector{Int64}
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
    ltmp_coral_cover::Array{Union{Missing, Float64}, 2}
    coral_composition::Array{Union{Missing, Float64}, 3}
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
function composition_idx_to_domain(loc_data_store::LocationDataStore, composition_idx::Int64)
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
include("./plotting/plotting.jl")
include("./plotting/progress.jl")
