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


global config = TOML.parsefile("calib_config.toml")
global path_configs = config["data_paths"]

global reefmod_domain_path = path_configs["reefmod_domain"]
global rme_domain_path = path_configs["rme_domain"]
global canonical_path = path_configs["canonical_path"]
global init_guess_fn = path_configs["init_guess_fn"]  # optional

global manta_tow_class_path = joinpath(pwd(), "../", path_configs["manta_tow_path"])  # target data for location classes
global historic_dhw_path = joinpath(rme_domain_path, "data_files", "dhw", "GBR_past_DHW_CRW_5km_1985_2022_Dec_2022.csv")
global ltmp_reef_data_path = joinpath(pwd(), "../", path_configs["ltmp_reef_data"])  # target data for ltmp locs
global composition_path = joinpath(pwd(), "../", path_configs["composition_netcdf"])

global ltmp_shp = joinpath(pwd(), "../", path_configs["ltmp_shp"])
global ltmp_modelled_obs = joinpath(pwd(), "../", path_configs["ltmp_modelled_obs"])
global classification_path = joinpath(pwd(), "../", path_configs["classification_path"])  # location classes
global init_cover_fn = joinpath(pwd(), "../", path_configs["init_cover_fn"])

global OUT_DIR = path_configs["out_dir"]

global start_year = 2008
global end_year = 2022
global n_scenarios = 16

if !@isdefined(OPTIONS)
    global uniform_initial_cover = false
    global rerun = false

    global reefmod_domain = true
    global reload_domain = false
    global reload_shp = false
    global reload_ltmp = false
    global reload_canonical = false

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
