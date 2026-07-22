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
    result_fn::String
    rng_seed::Int
end

"""
    load_config(config_path::String=<repo_root>/config.toml)::CalibrationConfig

Parse `config.toml` and resolve all calibration input/output paths, defaulting any path not
given explicitly to the corresponding file under `datasets/`.

Also sets `ENV["ADRIA_RNG_SEED"]` and `ENV["ADRIA_DEBUG"]` as a side effect: ADRIA's
`set_random_seed`/scenario setup read these env vars directly at `run_model` time with no
keyword-argument alternative, and `ADRIA.setup()`'s own `config.toml` resolution is
`pwd()`-relative and fails silently whenever the script isn't run from the repo root — so we
set `ADRIA_RNG_SEED` explicitly here rather than relying on it.
"""
function load_config(
    config_path::String=joinpath(dirname(dirname(@__DIR__)), "config.toml")
)::CalibrationConfig
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
    result_fn = get(output_config, "result_filename", "results.dat")

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
        result_fn,
        rng_seed
    )
end

"""
    load_domain(config::CalibrationConfig)

Load the RME domain and attach historical DHW/cyclone-mortality data used in place of
ADRIA's default projected scenarios.
"""
function load_domain(config::CalibrationConfig)
    @info "Loading RMEDomain"
    dom = ADRIA.load_domain(
        RMEDomain, config.rme_domain_path, "45", timeframe=(START_YEAR, END_YEAR)
    )

    @info "Attaching historic DHW and Cyclone/COTS data"
    new_dhw_scens = open_dataset(config.historical_dhw_path).dhw_scens
    dom.dhw_scens .= read(new_dhw_scens[timesteps=At(START_YEAR:END_YEAR), scenarios=1])

    new_cyclone_mortality_scens =
        open_dataset(config.historical_cyclone_mortality_path).disturbance_mortality_scens
    dom.cyclone_mortality_scens .= read(new_cyclone_mortality_scens[:, :, :, [1]])

    return dom
end

"""
    load_location_classification(loc_class_path::String)::DataFrame

Load the per-location consecutive-classification CSV consumed by `run_calibration`/
`construct_cover!`.
"""
function load_location_classification(loc_class_path::String)::DataFrame
    return CSV.read(loc_class_path, DataFrame)
end
