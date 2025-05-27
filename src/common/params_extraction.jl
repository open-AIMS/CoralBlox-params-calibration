# Coral params
# CORAL_PARAM_NAMES
calibrated_params[PARAM_IDXS[1]:PARAM_IDXS[2]]
coral_params_df = DataFrame(calibrated_params[PARAM_IDXS[1]:PARAM_IDXS[2]]', CORAL_PARAM_NAMES)
CSV.write(joinpath("C:/Users/pribeiro/AIMS", "coral_params.csv"), coral_params_df)

# Scale Factors params
scale_params_df = DataFrame(calibrated_params[PARAM_IDXS[3]:PARAM_IDXS[6]]', vcat(SCALE_FACTOR_NAMES, GROWTH_ACCEL_NAMES))
CSV.write(joinpath("C:/Users/pribeiro/AIMS", "scale_params.csv"), scale_params_df)
