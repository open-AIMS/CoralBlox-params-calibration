"""
    build_params_dataset(calibrated_params, param_idxs, coral_param_names, biogroups_ordering, growth_accel_names, depth_atten_names)::Dataset

Unpack a flat calibrated-parameter vector into a labelled YAXArray Dataset.

Reshape ordering follows the layout produced by `sc_fg_param_idxs` and
`accel_params_vec_to_array`:
- `linear_extension`, `mb_rate`: flat order is group-slowest, size-fastest → reshape
  `(n_sizes, n_groups)` then permute to `(n_groups, n_sizes)`
- `dist_mean`, `dist_std`: calibrated per functional group only, so each comes from its own
  parameter block as a length-`N_TAXA` vector rather than from the coral block
- `linear_extension_scale`, `mb_rate_scale`: flat order is group-slowest, biogroup-fastest
  → reshape `(n_biogroups, n_groups)` then permute to `(n_groups, n_biogroups)`
"""
function build_params_dataset(
    calibrated_params::Vector{Float64},
    param_idxs::Vector{Int64},
    coral_param_names,
    biogroups_ordering::Vector{Int64},
    growth_accel_names::Vector{String},
    depth_atten_names::Vector{String}
)::Dataset
    fgroup_names = ADRIA.functional_group_names()
    n_biogroups = length(biogroups_ordering)

    fg_axis = Dim{:functional_group}(string.(fgroup_names))
    sc_axis = Dim{:size_class}(1:N_SIZE_CLASSES)
    bg_axis = Dim{:cb_calib_group}(biogroups_ordering)
    ap_axis = Dim{:accel_param}(["steepness", "height", "midpoint"])
    # Labels must be ADRIA's DepthAttenuation field names - `_depth_attenuation_calib_overrides`
    # looks values up by this axis, not by position.
    da_axis = Dim{:depth_atten_param}(depth_atten_names)

    coral_params = calibrated_params[param_idxs[1]:param_idxs[2]]

    le_mask = occursin.(Ref("_linear_extension"), string.(coral_param_names))
    mb_mask = occursin.(Ref("_mb_rate"), string.(coral_param_names))
    les_mask = occursin.(Ref("linear_extension_scale"), string.(coral_param_names))
    mbs_mask = occursin.(Ref("mb_rate_scale"), string.(coral_param_names))

    _fg_sc(mask, desc, units=nothing) = YAXArray(
        (fg_axis, sc_axis),
        permutedims(reshape(coral_params[mask], N_SIZE_CLASSES, N_TAXA), (2, 1)),
        units === nothing ? Dict("description" => desc) : Dict("description" => desc, "units" => units)
    )

    # Per functional group, shared by every size class on the ADRIA side
    _fg(vals, desc) = YAXArray(
        (fg_axis,),
        vals,
        Dict(
            "description" => "$(desc), per functional group (shared by all size classes)",
            "units" => "DHW"
        )
    )

    _fg_bg(mask, desc) = YAXArray(
        (fg_axis, bg_axis),
        permutedims(reshape(coral_params[mask], n_biogroups, N_TAXA), (2, 1)),
        Dict("description" => desc)
    )

    accel_params = calibrated_params[param_idxs[3]:param_idxs[4]]
    steep = accel_params[occursin.(Ref("steepness"), growth_accel_names)]
    height = accel_params[occursin.(Ref("height"), growth_accel_names)]
    mid = accel_params[occursin.(Ref("midpoint"), growth_accel_names)]

    return Dataset(;
        linear_extension=_fg_sc(le_mask, "Linear extension rate", "m"),
        mb_rate=_fg_sc(mb_mask, "Background mortality rate"),
        dist_mean=_fg(
            calibrated_params[param_idxs[11]:param_idxs[12]],
            "DHW tolerance distribution mean"
        ),
        dist_std=_fg(
            calibrated_params[param_idxs[9]:param_idxs[10]],
            "DHW tolerance distribution standard deviation"
        ),
        linear_extension_scale=_fg_bg(les_mask, "Linear extension biogroup scale factor"),
        mb_rate_scale=_fg_bg(mbs_mask, "Background mortality biogroup scale factor"),
        growth_acceleration=YAXArray(
            (bg_axis, ap_axis),
            hcat(steep, height, mid),
            Dict("description" => "Logistic growth acceleration parameters per biogroup")
        ),
        depth_attenuation=YAXArray(
            (da_axis,),
            calibrated_params[param_idxs[7]:param_idxs[8]],
            Dict("description" => "GBR-wide depth attenuation parameters of surface DHW")
        ),
        properties=Dict(
            "description" => "Calibrated coral model parameters. Included in each variable (key)
            represents: linear extensions (linear_extension), background mortality rates
            (mb_rate), heat tolerance distribution means (dist_mean) and standard
            deviations (dist_std), linear extension scale factors per CB_GROUP
            (linear_extension_scale), background mortality rates scale factors per CB_GROUP
            (mb_rate_scale), growth acceleration steepness, height, and midpoints
            (growth_acceleration), and the GBR-wide depth attenuation of surface DHW
            (depth_attenuation)."
        )
    )
end

"""
    build_init_cover_dataset(dom, init_cover_sample, location_types, calibrated_params, param_idxs, observations, biogroup_ord)::Dataset

Compute the calibrated historical initial coral cover without mutating `dom`.

Applies `compute_cover` to get the GBR-wide historical cover from location-type samples,
then calls `insert_init_loc_cover!` on a temporary domain copy to override target-location
cover with photogrammetry-derived composition and calibrated size-class lambdas.
"""
function build_init_cover_dataset(
    dom,
    init_cover_sample::Vector{Float64},
    location_types::AbstractVector{Int64},
    calibrated_params::Vector{Float64},
    param_idxs::Vector{Int64},
    observations::LocationDataStore,
    biogroup_ord::Vector{Int64}
)::Dataset
    cover = compute_cover(init_cover_sample, location_types, ADRIA.n_locations(dom))
    tmp_dom = copy(dom)
    tmp_dom.init_coral_cover = copy(dom.init_coral_cover)
    tmp_dom.init_coral_cover .= cover
    insert_init_loc_cover!(tmp_dom, calibrated_params[param_idxs[5]:param_idxs[6]], observations, biogroup_ord)
    return Dataset(; tmp_dom.init_coral_cover)
end

"""
    export_calibration_products(dom, init_cover_sample, location_types, calibrated_params, param_idxs, coral_param_names, growth_accel_names, depth_atten_names, observations, biogroup_ord; out_dir)::String

Write the two NetCDF products of a calibration run to `{out_dir}/params/`:
`calibrated_params.nc` (ADRIA `calib_params_fn` input) and `historic_init_cover.nc`.

Returns the `params/` directory path. Existing files are overwritten.
"""
function export_calibration_products(
    dom,
    init_cover_sample::Vector{Float64},
    location_types::AbstractVector{Int64},
    calibrated_params::Vector{Float64},
    param_idxs::Vector{Int64},
    coral_param_names,
    growth_accel_names::Vector{String},
    depth_atten_names::Vector{String},
    observations::LocationDataStore,
    biogroup_ord::Vector{Int64};
    out_dir::String
)::String
    # ADRIA names its per-biogroup factors `*_cb_group_<i>_*` with `i` positional (1:n).
    # `load_calib_params` rebuilds those names from this file's `cb_calib_group` axis and
    # silently falls back to defaults on a miss, so a non-positional axis would produce an
    # uncalibrated domain with no error.
    @assert biogroup_ord == collect(1:length(biogroup_ord)) (
        "CB_CALIB_GROUPS must be 1:n to round-trip through ADRIA's calib_params_fn, got " *
        "$(biogroup_ord)"
    )

    params_dir_path = joinpath(out_dir, "params")
    mkpath(params_dir_path)

    savedataset(
        build_params_dataset(
            calibrated_params, param_idxs, coral_param_names, biogroup_ord,
            growth_accel_names, depth_atten_names
        );
        path=joinpath(params_dir_path, "calibrated_params.nc"),
        driver=:netcdf,
        overwrite=true
    )
    savedataset(
        build_init_cover_dataset(
            dom, init_cover_sample, location_types, calibrated_params, param_idxs,
            observations, biogroup_ord
        );
        path=joinpath(params_dir_path, "historic_init_cover.nc"),
        driver=:netcdf,
        overwrite=true
    )

    @info "Wrote calibration products to $(params_dir_path)"
    return params_dir_path
end
