# Result analysis and plot generation script.
#
# Prerequisites: run 01_a_setup.jl → 01_b_data_split.jl → 02_location_calibration.jl
# to produce the calibrated parameter file referenced by OUT_DIR / RESULT_FN.
#
# Using CoralBloxCalib first means LocationDataStore (and related helpers) come
# from the package, so the type is shared with the viz submodule without
# duplication.

using CoralBloxCalib
import CoralBloxCalib.viz

# ----- Package-level imports needed for setup ---------------------------------
using ADRIA
using ADRIA: GDF, AG
using CSV, DataFrames, YAXArrays, Serialization, TOML

# ----- Paths and config (mirrors common.jl without the global side-effects) ---

_script_dir = dirname(@__FILE__)
_src_dir    = joinpath(dirname(dirname(_script_dir)), "src")
_root_dir   = dirname(_src_dir)

CONFIG = TOML.parsefile(joinpath(_root_dir, "config.toml"))

DOMAIN_CONFIG       = CONFIG["calibration"]["domains"]
GEOSPATIAL_CONFIG   = CONFIG["calibration"]["geospatial"]
INITIALISATION_CONFIG = CONFIG["calibration"]["initialisation"]
OUTPUT_CONFIG       = CONFIG["calibration"]["outputs"]

RME_DOMAIN_PATH     = DOMAIN_CONFIG["rme_domain"]
CANONICAL_PATH      = GEOSPATIAL_CONFIG["canonical_path"]
LTMP_SHP_PATH       = get(GEOSPATIAL_CONFIG, "ltmp_shp",
    joinpath(_root_dir, "datasets", "spatial_data", "gbr_3Zone 2.shp"))
LOC_CLASS_PATH      = get(GEOSPATIAL_CONFIG, "classification_path",
    joinpath(_root_dir, "datasets", "spatial_data", "location_classification_MPA.csv"))

TARGET_CONFIG        = get(CONFIG["calibration"], "observations", Dict())
LTMP_REEF_DATA_PATH  = get(TARGET_CONFIG, "ltmp_reef_data",
    joinpath(_root_dir, "datasets", "ltmp_data", "manta_tow_data_reef_lvl.gpkg"))
COMPOSITION_PATH     = get(TARGET_CONFIG, "composition_netcdf",
    joinpath(_root_dir, "datasets", "ltmp_data", "coral_composition.nc"))
LTMP_MODELLED_OBS_PATH = get(TARGET_CONFIG, "ltmp_modelled_obs",
    joinpath(_root_dir, "datasets", "ltmp_data", "modelled_brms.beta.ry.disp.csv"))

INIT_COVER_PATH = get(INITIALISATION_CONFIG, "init_cover_filepath",
    joinpath(_root_dir, "datasets", "spatial_data", "init_cover.dat"))

HISTORICAL_DHW_PATH = joinpath(_src_dir, "historical_dhw_data_gen", "data", "historical_dhw.nc")
HISTORICAL_CYCLONE_MORTALITY_PATH = joinpath(
    _src_dir, "historical_disturbance_mortality_data_gen",
    "data", "historical_disturbance_mortality_rates",
    "historical_disturbance_mortality_rates.nc"
)

OUT_DIR   = OUTPUT_CONFIG["out_dir"]
RESULT_FN = get(OUTPUT_CONFIG, "result_filename", "results.dat")

ENV["ADRIA_RNG_SEED"] = string(CONFIG["operation"]["rng_seed"])

# ----- Setup: inlined from 03_a_analysis_setup.jl and 01_a_setup.jl ----------

using Statistics

@info "Loading RMEDomain"
dom = ADRIA.load_domain(RMEDomain, RME_DOMAIN_PATH, "45"; timeframe=(START_YEAR, END_YEAR))

@info "Attaching historic DHW and Cyclone/COTS data"
new_dhw_scens = open_dataset(HISTORICAL_DHW_PATH).dhw_scens
dom.dhw_scens .= read(new_dhw_scens[timesteps=At(START_YEAR:END_YEAR), scenarios=1])

new_cyc_scens = open_dataset(HISTORICAL_CYCLONE_MORTALITY_PATH).disturbance_mortality_scens
dom.cyclone_mortality_scens .= read(new_cyc_scens[:, :, :, [1]])

# Region masks
region_shps = GDF.read(LTMP_SHP_PATH)
function _region_shape_mask(dom, region_shapes, idx)::BitVector
    region_shape = region_shapes.geometry[idx]
    geoms = GDF.GeoInterface.geometry.(eachrow(dom.loc_data))
    return [AG.contains(region_shape, AG.centroid(geom)) for geom in geoms]
end
const NORTH_MASK   = _region_shape_mask(dom, region_shps, 1)
const CENTRAL_MASK = _region_shape_mask(dom, region_shps, 2)
const SOUTH_MASK   = _region_shape_mask(dom, region_shps, 3)

# LTMP modelled regional data
ltmp_data = CSV.read(LTMP_MODELLED_OBS_PATH, DataFrame, header=true)
ltmp_data[!, :Region] = String.(ltmp_data[:, :Region])

function _ltmp_region(region_name, df)
    mask = [region_name == r for r in df.Region]
    return df[mask, :]
end
ltmp_north, ltmp_central, ltmp_south = _ltmp_region.(
    ["Northern GBR", "Central GBR", "Southern GBR"], [ltmp_data]
)

location_classification = CSV.read(LOC_CLASS_PATH, DataFrame)

# ----- Build calibration / validation stores (mirrors 01_b_data_split.jl) ----

include(joinpath(_src_dir, "01_b_data_split.jl"))   # defines CALIBRATION_STORE, VALIDATION_STORE

# ----- 03_a: load cover, run model -------------------------------------------

init_cover = deserialize(INIT_COVER_PATH)
include(joinpath(_src_dir, "common", "cover_construction.jl"))
construct_cover!(dom, init_cover, location_classification.consecutive_classification)

calibrated_params = deserialize(joinpath(OUT_DIR, RESULT_FN))

include(joinpath(_src_dir, "common", "params_extraction.jl"))  # setup_run lives here if needed
# If setup_run is defined elsewhere, adjust the include path accordingly.
dom, scen = setup_run(dom, calibrated_params)
rs_raw = ADRIA.run_model(dom, scen[1, :])
mkpath(OUT_DIR)

# ----- Generate plots ---------------------------------------------------------

# Observation locations map
f_obs_loc_map = viz.plot_observation_locs(CALIBRATION_STORE, VALIDATION_STORE)
save(joinpath(OUT_DIR, "obs_loc_map.png"), f_obs_loc_map)

# Functional group cover proportions
f_taxa_cover = viz.taxa_cover_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_cov.png"), f_taxa_cover)

f_taxa_pop = viz.taxa_population_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_pop.png"), f_taxa_pop)

f_size_class = viz.temporal_size_class_proportions(rs_raw.raw; fig_size=(900, 600))
save(joinpath(OUT_DIR, "locs_size.png"), f_size_class)

# Summary stats
n_validation_locs = length(VALIDATION_STORE.ltmp_unique_ids)
validation_sccs = [
    viz.collect_error_stats(rs_raw.raw, id, dom; observations=VALIDATION_STORE).srcc
    for id in 1:n_validation_locs
]
validation_sccs_sortperm = sortperm(validation_sccs)
@info "Three highest SCC validation reefs: $(VALIDATION_STORE.ltmp_unique_ids[validation_sccs_sortperm][end-2:end])"
@info "Three lowest  SCC validation reefs: $(VALIDATION_STORE.ltmp_unique_ids[validation_sccs_sortperm][1:3])"

stats_calib = viz.collect_error_stats(rs_raw.raw, dom; observations=CALIBRATION_STORE)
stats_valid = viz.collect_error_stats(rs_raw.raw, dom; observations=VALIDATION_STORE)
n_calib = length(CALIBRATION_STORE.ltmp_unique_ids)
n_valid = length(VALIDATION_STORE.ltmp_unique_ids)
@info "Calibration locations where model outperforms benchmark: $(sum(stats_calib.rmse_model .< stats_calib.rmse_benchmark)) / $n_calib"
@info "Validation  locations where model outperforms benchmark: $(sum(stats_valid.rmse_model .< stats_valid.rmse_benchmark)) / $n_valid"
@info "Mean model calib. RMSE: $(mean(stats_calib.rmse_model))"
@info "Mean model valid.  RMSE: $(mean(stats_valid.rmse_model))"
@info "Mean model calib. SRCC: $(mean(stats_calib.srcc))"
@info "Mean model valid.  SRCC: $(mean(stats_valid.srcc))"

# Regional comparison plots (03_b)
viz.save_regional_analysis_plots(
    rs_raw.raw, dom, CALIBRATION_STORE, VALIDATION_STORE, OUT_DIR,
    ltmp_north, ltmp_central, ltmp_south,
    NORTH_MASK, CENTRAL_MASK, SOUTH_MASK
)

# Metric analysis plots (03_c)
viz.save_metric_analysis_plots(
    rs_raw.raw, dom, CALIBRATION_STORE, VALIDATION_STORE, OUT_DIR
)

# Per-location time-series plots (03_d)
cyc_scens = dom.cyclone_mortality_scens[scenarios=1]
dhw_scens  = dom.dhw_scens[scenarios=1]
disturbances = open_dataset(
    joinpath(_root_dir, "datasets", "ltmp_data", "disturbances.nc")
).layer

viz.save_location_timeseries_plots(
    rs_raw.raw, dom, CALIBRATION_STORE, VALIDATION_STORE, OUT_DIR,
    dhw_scens, cyc_scens, disturbances
)
