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

## Mb_rate scale factor
mb_rate_scale_mask = occursin.(Ref("mb_rate_scale"), string.(CORAL_PARAM_NAMES))
mb_rate_scale_data = calibrated_params[(PARAM_IDXS[1]:PARAM_IDXS[2])[mb_rate_scale_mask]]
mb_rate_scale_df = DataFrame(reshape(mb_rate_scale_data, 12, 5), fgroup_names)
CSV.write(joinpath(params_dir_path, "mb_rate_scale_params.csv"), mb_rate_scale_df)

# Growth acceleration params
steepness_mask = occursin.(Ref("steepness"), GROWTH_ACCEL_NAMES)
steepness_data = calibrated_params[PARAM_IDXS[3]:PARAM_IDXS[4]][steepness_mask]

height_mask = occursin.(Ref("height"), GROWTH_ACCEL_NAMES)
height_data = calibrated_params[PARAM_IDXS[3]:PARAM_IDXS[4]][height_mask]

midpoint_mask = occursin.(Ref("midpoint"), GROWTH_ACCEL_NAMES)
midpoint_data = calibrated_params[PARAM_IDXS[3]:PARAM_IDXS[4]][midpoint_mask]

growth_accel_params_df = DataFrame(hcat(steepness_data, height_data, midpoint_data), ["steepness", "height", "midpoint"])
CSV.write(joinpath(params_dir_path, "growth_accel_params.csv"), growth_accel_params_df)

# Historical initial cover
# TODO There's probably a better way of doing this but for now I'm just using what we already have
if !@isdefined(dom)
    include("../1_setup.jl")
    include("param_bounds.jl")
end
init_cover = deserialize(INIT_COVER_PATH)

construct_cover!(
    dom, init_cover, location_classification.consecutive_classification
)
coral_param_fn = joinpath(OUT_DIR, RESULT_FN)
calibrated_params = deserialize(coral_param_fn)

# Load values into scenario dataframe
dom, scen = setup_run(
    dom,
    calibrated_params
)

savedataset(
    Dataset(; dom.init_coral_cover);
    path=joinpath(params_dir_path, "historic_init_cover.nc"),
    driver=:netcdf,
    overwrite=true
)
