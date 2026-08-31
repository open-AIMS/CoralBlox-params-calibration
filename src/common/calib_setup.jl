"""
    CalibrationConfig

Explicit config/path bundle produced by [`load_config`](@ref), replacing the ~20 `global`
constants previously populated by the old `src/common/setup.jl` include-time script.
"""
struct CalibrationConfig
    rme_domain_path::String
    historical_dhw_path::String
    historical_cyclone_mortality_path::String
    canonical_path::String
    ltmp_shp_path::String
    loc_class_path::String
    ltmp_reef_data_path::String
    composition_path::String
    ltmp_modelled_obs_path::String
    init_cover_path::String
    init_guess_path::String
    out_dir::String
    rng_seed::Int
end

"""
    load_config(config_path::String)::CalibrationConfig

Parse the `config.toml` at `config_path` and resolve all calibration input/output paths,
defaulting any path not given explicitly to the corresponding file under this package's
`datasets/` directory.

`config_path` has no default: it must be supplied by the caller (typically
`joinpath(@__DIR__, "config.toml")` resolved against the *caller's* own location) rather than
resolved relative to this package's install location, so `load_config` behaves correctly when
called from outside this repo.

Also sets `ENV["ADRIA_RNG_SEED"]` and `ENV["ADRIA_DEBUG"]` as a side effect: ADRIA's
`set_random_seed`/scenario setup read these env vars directly at `run_model` time with no
keyword-argument alternative, and `ADRIA.setup()`'s own `config.toml` resolution is
`pwd()`-relative and fails silently whenever the script isn't run from the repo root — so we
set `ADRIA_RNG_SEED` explicitly here rather than relying on it.
"""
function load_config(config_path::String)::CalibrationConfig
    src_path = dirname(@__DIR__)
    root_path = dirname(src_path)
    datasets_path = joinpath(root_path, "datasets")

    config = TOML.parsefile(config_path)

    domain_config = config["calibration"]["domains"]
    geospatial_config = config["calibration"]["geospatial"]
    initialisation_config = get(config["calibration"], "initialisation", Dict())
    output_config = config["calibration"]["outputs"]
    target_config = get(config["calibration"], "observations", Dict())

    # Single source of truth for the run seed, shared by ADRIA and BlackBoxOptim (passed
    # explicitly as RngSeed).
    rng_seed = config["operation"]["rng_seed"]
    ENV["ADRIA_RNG_SEED"] = string(rng_seed)
    ENV["ADRIA_DEBUG"] = false

    rme_domain_path = domain_config["rme_domain"]
    historical_dhw_path = joinpath(
        src_path, "historical_dhw_data_gen", "data", "historical_dhw.nc"
    )
    historical_cyclone_mortality_path = joinpath(
        src_path,
        "historical_disturbance_mortality_data_gen/",
        "data/",
        "historical_disturbance_mortality_rates/",
        "historical_disturbance_mortality_rates.nc"
    )

    canonical_path = geospatial_config["canonical_path"]
    ltmp_shp_path = get(
        geospatial_config, "ltmp_shp", joinpath(datasets_path, "spatial_data/gbr_3Zone 2.shp")
    )
    loc_class_path = get(
        geospatial_config,
        "classification_path",
        joinpath(datasets_path, "spatial_data/location_classification_MPA.csv")
    )

    ltmp_reef_data_path = get(
        target_config,
        "ltmp_reef_data",
        joinpath(datasets_path, "ltmp_data/manta_tow_data_reef_lvl.gpkg")
    )
    composition_path = get(
        target_config,
        "composition_netcdf",
        joinpath(datasets_path, "ltmp_data/coral_composition.nc")
    )
    ltmp_modelled_obs_path = get(
        target_config,
        "ltmp_modelled_obs",
        joinpath(datasets_path, "ltmp_data/modelled_brms.beta.ry.disp.csv")
    )

    init_cover_path = get(
        initialisation_config,
        "init_cover_filepath",
        joinpath(datasets_path, "spatial_data/init_cover.dat")
    )
    init_guess_path = get(initialisation_config, "init_guess_filepath", "")

    out_dir = output_config["out_dir"]

    return CalibrationConfig(
        rme_domain_path,
        historical_dhw_path,
        historical_cyclone_mortality_path,
        canonical_path,
        ltmp_shp_path,
        loc_class_path,
        ltmp_reef_data_path,
        composition_path,
        ltmp_modelled_obs_path,
        init_cover_path,
        init_guess_path,
        out_dir,
        rng_seed
    )
end

"""
    load_domain(config::CalibrationConfig; calib_params_fn::String="")

Load the RME domain and attach historical DHW/cyclone-mortality data used in place of
ADRIA's default projected scenarios.

`calib_params_fn` is passed through to ADRIA: given a `calibrated_params.nc` written by
`export_calibration_products`, the returned domain's model spec already carries the
calibrated coral, growth-acceleration and depth-attenuation values, so
`ADRIA.param_table(dom)` yields a calibrated scenario without going through `setup_run`.
"""
function load_domain(config::CalibrationConfig; calib_params_fn::String="")
    @info "Loading RMEDomain"
    dom = ADRIA.load_domain(
        RMEDomain, config.rme_domain_path, "45", timeframe=(START_YEAR, END_YEAR),
        calib_params_fn=calib_params_fn
    )

    @info "Attaching historic DHW and Cyclone/COTS data"
    new_dhw_scens = open_dataset(config.historical_dhw_path).dhw_scens
    dom.dhw_scens .= read(new_dhw_scens[timesteps=At(START_YEAR:END_YEAR), scenarios=1])

    new_cyclone_mortality_scens =
        open_dataset(config.historical_cyclone_mortality_path).disturbance_mortality_scens
    dom.cyclone_mortality_scens .= read(new_cyclone_mortality_scens[:, :, :, [1]])

    @info "Overwriting location classification/identity data from canonical geopackage"
    _attach_canonical_loc_data!(dom, config.canonical_path)

    return dom
end

"""
    _attach_canonical_loc_data!(dom::RMEDomain, canonical_path::String)::Nothing

Overwrite `dom.loc_data`'s reef classification/identity columns (`UNIQUE_ID`, `GBRMPA_ID`,
`management_area_short`, `GBRMPA_BIOREGION`, `CB_CALIB_GROUPS`, `k`, `area`, `cluster_id`)
with values read fresh from the canonical geopackage at `canonical_path`, keyed by
`dom.loc_ids`.

These columns are otherwise baked into the RME domain package's own `reefmod_gbr.gpkg` at
whatever time that package was built, giving two independently-updatable copies of the same
source-of-truth reef data. Sourcing them from `canonical_path` on every load removes that
drift risk. Physical/simulation-only columns (`depth_med`, `mean_to_neighbor`, `zone_type`)
are untouched: they either have no canonical equivalent or are deliberately RME-specific
(e.g. `depth_med` is flattened to a constant 7.0 by ADRIA's RME domain loader).

Fails loudly (rather than silently producing `missing`/`NaN` columns) if `canonical_path`
is missing a required source column or a location id present in `dom.loc_ids`.
"""
function _attach_canonical_loc_data!(dom::RMEDomain, canonical_path::String)::Nothing
    canonical_gdf = GDF.read(canonical_path)

    required_cols = (
        dom.loc_id_col, "UNIQUE_ID", "GBRMPA_ID", "management_area_short",
        "GBRMPA_BIOREGION", "CB_CALIB_GROUPS", "ReefMod_habitable_proportion",
        "ReefMod_area_m2", "reef_name"
    )
    missing_cols = setdiff(required_cols, names(canonical_gdf))
    if !isempty(missing_cols)
        error(
            "canonical geopackage at $(canonical_path) is missing required column(s): " *
            "$(join(missing_cols, ", "))"
        )
    end

    # `dom.loc_ids` (from the domain's DHW dataset) is keyed by `dom.loc_id_col`
    # (e.g. "RME_UNIQUE_ID"), which is NOT always identical to the geopackage's own
    # `UNIQUE_ID` column row-for-row - join on the id column the domain actually uses.
    join_col = canonical_gdf[!, dom.loc_id_col]
    row_idx = [findfirst(==(id), join_col) for id in dom.loc_ids]
    missing_mask = isnothing.(row_idx)
    if any(missing_mask)
        missing_ids = dom.loc_ids[missing_mask]
        error(
            "canonical geopackage at $(canonical_path) is missing " *
            "$(length(missing_ids)) location id(s) present in the RME domain, " *
            "e.g. $(first(missing_ids, min(5, length(missing_ids))))"
        )
    end
    canon = canonical_gdf[row_idx, :]

    dom.loc_data[:, :UNIQUE_ID] = canon.UNIQUE_ID
    dom.loc_data[:, :GBRMPA_ID] = canon.GBRMPA_ID
    dom.loc_data[:, :management_area_short] = canon.management_area_short
    dom.loc_data[:, :GBRMPA_BIOREGION] = canon.GBRMPA_BIOREGION
    dom.loc_data[:, :CB_CALIB_GROUPS] = canon.CB_CALIB_GROUPS
    dom.loc_data[:, :k] = canon.ReefMod_habitable_proportion
    dom.loc_data[:, :area] = canon.ReefMod_area_m2
    dom.loc_data[:, :cluster_id] = canon.reef_name

    return nothing
end

"""
    load_location_classification(loc_class_path::String)::DataFrame

Load the per-location consecutive-classification CSV consumed by `run_calibration`/
`construct_cover!`.
"""
function load_location_classification(loc_class_path::String)::DataFrame
    return CSV.read(loc_class_path, DataFrame)
end
