using YAXArrays

include("param_bounds.jl")

coral_param_fn = joinpath(OUT_DIR, RESULT_FN)
calibrated_params = deserialize(coral_param_fn)

fgroup_names = ADRIA.functional_group_names()

# Dimension axes — match ADRIA array shapes from coral_factors.jl / GrowthAcceleration.jl
fg_axis = Dim{:functional_group}(string.(fgroup_names))
sc_axis = Dim{:size_class}(1:n_size_classes)
bg_axis = Dim{:cb_calib_group}(BIOGROUPS_ORDERING)
ap_axis = Dim{:accel_param}(["steepness", "height", "midpoint"])

coral_params_range = PARAM_IDXS[1]:PARAM_IDXS[2]

# Coral params — shapes match ADRIA: (n_groups=5, n_sizes=7)
# Data from sc_fg_param_idxs is ordered group-slowest, size-fastest → reshape (n_sizes, n_groups) then permute

# Linear extension
linear_extension_mask = occursin.(Ref("_linear_extension"), string.(CORAL_PARAM_NAMES))
linear_extension_data = calibrated_params[coral_params_range[linear_extension_mask]]
linear_extension_yax = YAXArray(
    (fg_axis, sc_axis),
    permutedims(reshape(linear_extension_data, n_size_classes, n_groups), (2, 1)),
    Dict("description" => "Linear extension rate", "units" => "m")
)

# Mb_rate
mb_rate_mask = occursin.(Ref("_mb_rate"), string.(CORAL_PARAM_NAMES))
mb_rate_data = calibrated_params[coral_params_range[mb_rate_mask]]
mb_rate_yax = YAXArray(
    (fg_axis, sc_axis),
    permutedims(reshape(mb_rate_data, n_size_classes, n_groups), (2, 1)),
    Dict("description" => "Background mortality rate")
)

# Dist mean
# dist_mean params are ordered size-slowest, group-fastest → reshape (n_groups, n_sizes) directly
dist_mean_mask = occursin.(Ref("_dist_mean"), string.(CORAL_PARAM_NAMES))
dist_mean_data = calibrated_params[coral_params_range[dist_mean_mask]]
dist_mean_yax = YAXArray(
    (fg_axis, sc_axis),
    reshape(dist_mean_data, n_groups, n_size_classes),
    Dict("description" => "DHW tolerance distribution mean")
)

# Linear extension scale factor — shape matches ADRIA: (n_groups=5, n_cb_calib_groups=12)
# Data ordered group-fastest, biogroup-slowest → reshape (n_biogroups, n_groups) then permute
linear_extension_scale_mask = occursin.(Ref("linear_extension_scale"), string.(CORAL_PARAM_NAMES))
linear_extension_scale_data = calibrated_params[coral_params_range[linear_extension_scale_mask]]
linear_extension_scale_yax = YAXArray(
    (fg_axis, bg_axis),
    permutedims(reshape(linear_extension_scale_data, n_biogroups, n_groups), (2, 1)),
    Dict("description" => "Linear extension biogroup scale factor")
)

# Mb_rate scale factor
mb_rate_scale_mask = occursin.(Ref("mb_rate_scale"), string.(CORAL_PARAM_NAMES))
mb_rate_scale_data = calibrated_params[coral_params_range[mb_rate_scale_mask]]
mb_rate_scale_yax = YAXArray(
    (fg_axis, bg_axis),
    permutedims(reshape(mb_rate_scale_data, n_biogroups, n_groups), (2, 1)),
    Dict("description" => "Background mortality biogroup scale factor")
)

# Growth acceleration params — shape matches ADRIA accel_params_vec_to_array: (n_biogroups=12, 3)
steepness_mask = occursin.(Ref("steepness"), GROWTH_ACCEL_NAMES)
height_mask = occursin.(Ref("height"), GROWTH_ACCEL_NAMES)
midpoint_mask = occursin.(Ref("midpoint"), GROWTH_ACCEL_NAMES)

growth_accel_range = PARAM_IDXS[3]:PARAM_IDXS[4]
steepness_data = calibrated_params[growth_accel_range][steepness_mask]
height_data = calibrated_params[growth_accel_range][height_mask]
midpoint_data = calibrated_params[growth_accel_range][midpoint_mask]

growth_accel_yax = YAXArray(
    (bg_axis, ap_axis),
    hcat(steepness_data, height_data, midpoint_data),
    Dict("description" => "Logistic growth acceleration parameters per biogroup")
)

# Save all parameters to a single NetCDF
ds = Dataset(;
    linear_extension=linear_extension_yax,
    mb_rate=mb_rate_yax,
    dist_mean=dist_mean_yax,
    linear_extension_scale=linear_extension_scale_yax,
    mb_rate_scale=mb_rate_scale_yax,
    growth_acceleration=growth_accel_yax,
    properties=Dict(
        "description" => "Calibrated coral model parameters. Included in each variable (key)
        represents: linear extensions (linear_extension), background mortality rates
        (mb_rate_scale), heat tolerance distribution means (dist_mean), linear extension
        scale factors per CB_GROUP (linear_extension_scale), background mortality rates scale
        factors per CB_GROUP (mb_rate_scale), growth acceleration steepness, height,
        and midpoints (growth_acceleration)."
    )
)

params_dir_path = joinpath(OUT_DIR, "params")
!isdir(params_dir_path) && mkdir(params_dir_path)

savedataset(
    ds;
    path=joinpath(params_dir_path, "calibrated_params.nc"),
    driver=:netcdf,
    overwrite=true
)

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
