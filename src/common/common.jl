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


global CONFIG = TOML.parsefile("calib_config.toml")

# Configuration sections
global DOMAIN_CONFIG         = CONFIG["Domains"]
global GEOSPATIAL_CONFIG     = CONFIG["Geospatial"]
global TARGET_CONFIG         = CONFIG["Target"]
global INITIALISATION_CONFIG = CONFIG["Initialisation"]
global OUTPUT_CONFIG         = CONFIG["Outputs"]

# ADRIA Domain paths
global REEFMOD_DOMAIN_PATH = DOMAIN_CONFIG["reefmod_domain"]
global RME_DOMAIN_PATH     = DOMAIN_CONFIG["rme_domain"]
global HISTORIC_DHW_PATH   = joinpath(RME_DOMAIN_PATH, "data_files", "dhw", "GBR_past_DHW_CRW_5km_1985_2022_Dec_2022.csv")

# Geospatial filepaths
global CANONICAL_PATH = GEOSPATIAL_CONFIG["canonical_path"]
global LTMP_SHP_PATH  = GEOSPATIAL_CONFIG["ltmp_shp"]
global LOC_CLASS_PATH = GEOSPATIAL_CONFIG["classification_path"] # location classes

# Calibration Target / Observational Data
global LOC_CLASS_TARGET_PATH  = TARGET_CONFIG["manta_tow_path"]  # target data for location classes
global LTMP_REEF_DATA_PATH    = TARGET_CONFIG["ltmp_reef_data"]  # target data for ltmp locs
global COMPOSITION_PATH       = TARGET_CONFIG["composition_netcdf"]
global LTMP_MODELLED_OBS_PATH = TARGET_CONFIG["ltmp_modelled_obs"]

# Initialisation filepaths
global INIT_COVER_PATH  = INITIALISATION_CONFIG["init_cover_filepath"]
global INIT_GUESS_PATH  = INITIALISATION_CONFIG["init_guess_filepath"]  # optiona;

# Output filepaths
global OUT_DIR   = OUTPUT_CONFIG["out_dir"]
global RESULT_FN = OUTPUT_CONFIG["result_filename"]

global start_year = 2008
global end_year = 2022
global n_scenarios = 16

if !@isdefined(OPTIONS)
    global reload_domain = false

    # define OPTIONS to prevent reinitialise on every include
    global OPTIONS = true
end

include("./perf_metrics.jl")
include("./plotting/plotting.jl")
include("./plotting/progress.jl")


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
