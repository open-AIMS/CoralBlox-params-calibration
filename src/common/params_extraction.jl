"""
    build_params_dataset(calibrated_params, param_idxs, coral_param_names, biogroups_ordering, growth_accel_names)::Dataset

Unpack a flat calibrated-parameter vector into a labelled YAXArray Dataset.

Reshape ordering follows the layout produced by `sc_fg_param_idxs` and
`accel_params_vec_to_array`:
- `linear_extension`, `mb_rate`, `dist_mean`, `dist_std`: flat order is group-slowest,
  size-fastest → reshape `(n_sizes, n_groups)` then permute to `(n_groups, n_sizes)`
- `linear_extension_scale`, `mb_rate_scale`: flat order is group-slowest, biogroup-fastest
  → reshape `(n_biogroups, n_groups)` then permute to `(n_groups, n_biogroups)`
"""
function build_params_dataset(
    calibrated_params::Vector{Float64},
    param_idxs::Vector{Int64},
    coral_param_names,
    biogroups_ordering::Vector{Int64},
    growth_accel_names::Vector{String}
)::Dataset
    fgroup_names = ADRIA.functional_group_names()
    n_biogroups = length(biogroups_ordering)

    fg_axis = Dim{:functional_group}(string.(fgroup_names))
    sc_axis = Dim{:size_class}(1:N_SIZE_CLASSES)
    bg_axis = Dim{:cb_calib_group}(biogroups_ordering)
    ap_axis = Dim{:accel_param}(["steepness", "height", "midpoint"])

    coral_params = calibrated_params[param_idxs[1]:param_idxs[2]]

    le_mask = occursin.(Ref("_linear_extension"), string.(coral_param_names))
    mb_mask = occursin.(Ref("_mb_rate"), string.(coral_param_names))
    dm_mask = occursin.(Ref("_dist_mean"), string.(coral_param_names))
    ds_mask = occursin.(Ref("_dist_std"), string.(coral_param_names))
    les_mask = occursin.(Ref("linear_extension_scale"), string.(coral_param_names))
    mbs_mask = occursin.(Ref("mb_rate_scale"), string.(coral_param_names))

    _fg_sc(mask, desc, units=nothing) = YAXArray(
        (fg_axis, sc_axis),
        permutedims(reshape(coral_params[mask], N_SIZE_CLASSES, N_TAXA), (2, 1)),
        units === nothing ? Dict("description" => desc) : Dict("description" => desc, "units" => units)
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
        dist_mean=_fg_sc(dm_mask, "DHW tolerance distribution mean", "DHW"),
        dist_std=_fg_sc(ds_mask, "DHW tolerance distribution standard deviation", "DHW"),
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
            (mb_rate), heat tolerance distribution means (dist_mean) and standard
            deviations (dist_std), linear extension scale factors per CB_GROUP
            (linear_extension_scale), background mortality rates scale factors per CB_GROUP
            (mb_rate_scale), growth acceleration steepness, height, and midpoints
            (growth_acceleration)."
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
