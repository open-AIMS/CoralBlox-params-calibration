using Revise, Infiltrator
using Serialization

using TOML
using Dates
using CSV, DataFrames, YAXArrays
using StatsBase, Statistics

using ADRIA
using ADRIA: GDF, AG, DimensionalData

src_path = dirname(@__DIR__)
root_path = dirname(src_path)
datasets_path = joinpath(dirname(src_path), "datasets")

global CONFIG = TOML.parsefile(joinpath(dirname(src_path), "config.toml"))

# Configuration sections
global DOMAIN_CONFIG = CONFIG["calibration"]["domains"]
global GEOSPATIAL_CONFIG = CONFIG["calibration"]["geospatial"]
global INITIALISATION_CONFIG = get(CONFIG["calibration"], "initialisation", Dict())
global OUTPUT_CONFIG = CONFIG["calibration"]["outputs"]

# Single source of truth for the run seed, shared by ADRIA and BlackBoxOptim (passed
# explicitly as RngSeed). We set ENV["ADRIA_RNG_SEED"] ourselves here rather than relying
# on ADRIA.setup() to wire it up: ADRIA.setup() resolves config.toml via
# joinpath(pwd(), "config.toml"), which fails silently (falls back to defaults) whenever
# pwd() isn't the repo root — e.g. when scripts are run from inside src/, per this repo's
# README. Setting the env var explicitly here makes ADRIA's RNG seeding correct
# regardless of pwd() or whether ADRIA.setup()'s own file lookup succeeds.
global RNG_SEED = CONFIG["operation"]["rng_seed"]
ENV["ADRIA_RNG_SEED"] = string(RNG_SEED)

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
global TARGET_CONFIG = get(CONFIG["calibration"], "observations", Dict())

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

if !@isdefined(OPTIONS)
    global reload_domain = false

    # define OPTIONS to prevent reinitialise on every include
    global OPTIONS = true
end


