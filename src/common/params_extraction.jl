include("param_bounds.jl")

coral_param_fn = joinpath(OUT_DIR, RESULT_FN)
calibrated_params = deserialize(coral_param_fn)

params_dir_path = joinpath(OUT_DIR, "params")
!isdir(params_dir_path) && mkdir(params_dir_path)

fgroup_names = ADRIA.functional_group_names()

# Coral params

## Linear extension
linear_extension_mask = occursin.(Ref("_linear_extension"), string.(CORAL_PARAM_NAMES))
linear_extension_data = calibrated_params[(PARAM_IDXS[1]:PARAM_IDXS[2])[linear_extension_mask]]
linear_extension_df = DataFrame(reshape(linear_extension_data, 7, 5), fgroup_names)
CSV.write(joinpath(params_dir_path, "linear_extension_params.csv"), linear_extension_df)

## Mb_rate
mb_rate_mask = occursin.(Ref("_mb_rate"), string.(CORAL_PARAM_NAMES))
mb_rate_data = calibrated_params[(PARAM_IDXS[1]:PARAM_IDXS[2])[mb_rate_mask]]
mb_rate_df = DataFrame(reshape(mb_rate_data, 7, 5), fgroup_names)
CSV.write(joinpath(params_dir_path, "mb_rate_params.csv"), mb_rate_df)

## Dist mean
dist_mean_mask = occursin.(Ref("_dist_mean"), string.(CORAL_PARAM_NAMES))
dist_mean_data = calibrated_params[(PARAM_IDXS[1]:PARAM_IDXS[2])[dist_mean_mask]]
dist_mean_df = DataFrame(reshape(dist_mean_data, 7, 5)', fgroup_names)
CSV.write(joinpath(params_dir_path, "dist_mean_params.csv"), dist_mean_df)

## Linear extension scale factor
linear_extension_scale_mask = occursin.(Ref("linear_extension_scale"), string.(CORAL_PARAM_NAMES))
linear_extension_scale_data = calibrated_params[(PARAM_IDXS[1]:PARAM_IDXS[2])[linear_extension_scale_mask]]
linear_extension_scale_df = DataFrame(reshape(linear_extension_scale_data, 12, 5), fgroup_names)
CSV.write(joinpath(params_dir_path, "linear_extension_scale_params.csv"), linear_extension_scale_df)

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
