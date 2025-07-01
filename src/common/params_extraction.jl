include("param_bounds.jl")

coral_param_fn = joinpath(OUT_DIR, RESULT_FN)
calibrated_params = deserialize(coral_param_fn)

params_dir_path = joinpath(OUT_DIR, "params")
!isdir(params_dir_path) && mkdir(params_dir_path)

# Coral params
# CORAL_PARAM_NAMES
coral_params_df = DataFrame(calibrated_params[PARAM_IDXS[1]:PARAM_IDXS[2]]', CORAL_PARAM_NAMES)
CSV.write(joinpath(params_dir_path, "coral_params.csv"), coral_params_df)

# Scale Factors params
scale_params_df = DataFrame(calibrated_params[PARAM_IDXS[3]:PARAM_IDXS[4]]', SCALE_FACTOR_NAMES)
CSV.write(joinpath(params_dir_path, "scale_params.csv"), scale_params_df)

growth_accel_params_df = DataFrame(calibrated_params[PARAM_IDXS[5]:PARAM_IDXS[6]]', GROWTH_ACCEL_NAMES)
CSV.write(joinpath(params_dir_path, "growth_accel_params.csv"), growth_accel_params_df)

# In case we want to check the calibrated values
# scale_lin_ext = calibrated_params[PARAM_IDXS[3]:PARAM_IDXS[4]][[occursin(r"linear_extension", s) for s in SCALE_FACTOR_NAMES]]
# mean(scale_lin_ext)
# std(scale_lin_ext)

# scale_mb_rate = calibrated_params[PARAM_IDXS[3]:PARAM_IDXS[4]][[occursin(r"mb_rate", s) for s in SCALE_FACTOR_NAMES]]
# mean(scale_mb_rate)
# std(scale_mb_rate)
