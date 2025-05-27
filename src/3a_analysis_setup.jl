init_cover = deserialize(INIT_COVER_PATH)

construct_cover!(
    dom, init_cover, location_classification.consecutive_classification
)

# ----- LOAD CALIBRATED RESULTS -----

coral_param_fn = joinpath(OUT_DIR, RESULT_FN)
calibrated_params = deserialize(coral_param_fn)

# Load values into scenario dataframe
dom, scen = setup_run(
    dom,
    calibrated_params
)

rs_raw = ADRIA.run_model(dom, scen[1, :])
