using YAXArrays

include("param_bounds.jl")

"""
    build_params_dataset(calibrated_params)::Dataset

Unpack a flat calibrated-parameter vector into a labelled YAXArray Dataset.

Reshape ordering follows the layout produced by `sc_fg_param_idxs` and
`ADRIA.accel_params_vec_to_array`:
- `linear_extension`, `mb_rate`: flat order is group-slowest, size-fastest
  → reshape `(n_sizes, n_groups)` then permute to `(n_groups, n_sizes)`
- `dist_mean`: flat order is size-slowest, group-fastest
  → reshape `(n_groups, n_sizes)` directly
- `linear_extension_scale`, `mb_rate_scale`: flat order is group-fastest, biogroup-slowest
  → reshape `(n_biogroups, n_groups)` then permute to `(n_groups, n_biogroups)`
"""
function build_params_dataset(calibrated_params::Vector{Float64})::Dataset
    fgroup_names = ADRIA.functional_group_names()

    fg_axis = Dim{:functional_group}(string.(fgroup_names))
    sc_axis = Dim{:size_class}(1:n_size_classes)
    bg_axis = Dim{:cb_calib_group}(BIOGROUPS_ORDERING)
    ap_axis = Dim{:accel_param}(["steepness", "height", "midpoint"])

    coral_params = calibrated_params[PARAM_IDXS[1]:PARAM_IDXS[2]]

    le_mask = occursin.(Ref("_linear_extension"), string.(CORAL_PARAM_NAMES))
    mb_mask = occursin.(Ref("_mb_rate"), string.(CORAL_PARAM_NAMES))
    dm_mask = occursin.(Ref("_dist_mean"), string.(CORAL_PARAM_NAMES))
    les_mask = occursin.(Ref("linear_extension_scale"), string.(CORAL_PARAM_NAMES))
    mbs_mask = occursin.(Ref("mb_rate_scale"), string.(CORAL_PARAM_NAMES))

    _fg_sc(mask, desc, units=nothing) = YAXArray(
        (fg_axis, sc_axis),
        permutedims(reshape(coral_params[mask], n_size_classes, n_groups), (2, 1)),
        units === nothing ? Dict("description" => desc) : Dict("description" => desc, "units" => units)
    )

    _fg_bg(mask, desc) = YAXArray(
        (fg_axis, bg_axis),
        permutedims(reshape(coral_params[mask], n_biogroups, n_groups), (2, 1)),
        Dict("description" => desc)
    )

    accel_params = calibrated_params[PARAM_IDXS[3]:PARAM_IDXS[4]]
    steep = accel_params[occursin.(Ref("steepness"), GROWTH_ACCEL_NAMES)]
    height = accel_params[occursin.(Ref("height"), GROWTH_ACCEL_NAMES)]
    mid = accel_params[occursin.(Ref("midpoint"), GROWTH_ACCEL_NAMES)]

    return Dataset(;
        linear_extension=_fg_sc(le_mask, "Linear extension rate", "m"),
        mb_rate=_fg_sc(mb_mask, "Background mortality rate"),
        dist_mean=YAXArray(
            (fg_axis, sc_axis),
            reshape(coral_params[dm_mask], n_groups, n_size_classes),
            Dict("description" => "DHW tolerance distribution mean")
        ),
        linear_extension_scale=_fg_bg(les_mask, "Linear extension biogroup scale factor"),
        mb_rate_scale=_fg_bg(mbs_mask, "Background mortality biogroup scale factor"),
        growth_acceleration=YAXArray(
            (bg_axis, ap_axis),
            hcat(steep, height, mid),
            Dict("description" => "Logistic growth acceleration parameters per biogroup")
        ),
        properties=Dict(
            "description" => "Calibrated coral model parameters. Included in each variable (key)
            represents: linear extensions (linear_extension), background mortality rates
            (mb_rate_scale), heat tolerance distribution means (dist_mean), linear extension
            scale factors per CB_GROUP (linear_extension_scale), background mortality rates scale
            factors per CB_GROUP (mb_rate_scale), growth acceleration steepness, height,
            and midpoints (growth_acceleration)."
        )
    )
end

calibrated_params = deserialize(joinpath(OUT_DIR, RESULT_FN))

params_dir_path = joinpath(OUT_DIR, "params")
!isdir(params_dir_path) && mkdir(params_dir_path)

savedataset(
    build_params_dataset(calibrated_params);
    path=joinpath(params_dir_path, "calibrated_params.nc"),
    driver=:netcdf,
    overwrite=true
)

"""
    build_init_cover_dataset(dom, init_cover_sample, location_types, calibrated_params)::Dataset

Compute the calibrated historical initial coral cover without mutating `dom`.

Applies `compute_cover` to get the GBR-wide historical cover from location-type samples,
then calls `insert_init_loc_cover!` on a temporary domain copy to override target-location
cover with photogrammetry-derived composition and calibrated size-class lambdas.
"""
function build_init_cover_dataset(
    dom,
    init_cover_sample::Vector{Float64},
    location_types::AbstractVector{Int64},
    calibrated_params::Vector{Float64}
)::Dataset
    cover = compute_cover(init_cover_sample, location_types, ADRIA.n_locations(dom))
    tmp_dom = copy(dom)
    tmp_dom.init_coral_cover = copy(dom.init_coral_cover)
    tmp_dom.init_coral_cover .= cover
    insert_init_loc_cover!(tmp_dom, calibrated_params[PARAM_IDXS[5]:PARAM_IDXS[6]])
    return Dataset(; tmp_dom.init_coral_cover)
end

if !@isdefined(dom)
    include("../1_setup.jl")
    include("param_bounds.jl")
end
init_cover_sample = deserialize(INIT_COVER_PATH)
calibrated_params = deserialize(joinpath(OUT_DIR, RESULT_FN))

savedataset(
    build_init_cover_dataset(dom, init_cover_sample, location_classification.consecutive_classification, calibrated_params);
    path=joinpath(params_dir_path, "historic_init_cover.nc"),
    driver=:netcdf,
    overwrite=true
)
