struct CalibConfig
    sample_bounds::Vector{Tuple{Float64,Float64}}
    biogroups_ordering::Vector{Int64}
    param_idxs::Vector{Int64}
    coral_param_names::Vector{Symbol}
    growth_accel_names::Vector{String}
end

function CalibConfig(dom)::CalibConfig
    linear_extension_mean = [
        0.0245124 0.0507098 0.0501524 0.0637291 0.0672375 0.0779938 0.0
        0.024296  0.0286086 0.0278468 0.0285766 0.0314185 0.0364447 0.0
        0.0175738 0.0168217 0.0150495 0.0165701 0.0165701 0.0165701 0.0
        0.0116048 0.00747208 0.00748131 0.00942616 0.0133995 0.0141176 0.0
        0.011934  0.00747208 0.00748131 0.00942616 0.0133995 0.0141176 0.0
    ]
    lin_ext_lb = flatten_group_size(linear_extension_mean[:, :] .* LIN_EXT_LB_FACTOR)
    lin_ext_ub = flatten_group_size(linear_extension_mean[:, :] .* LIN_EXT_UB_FACTOR)

    mb_rate_mean = [
        0.312661  0.194444  0.211039   0.192857   0.157895   0.142857   0.142857
        0.223847  0.130748  0.0915385  0.123348   0.110294   0.110294   0.110294
        0.218824  0.128571  0.078534   0.0833333  0.0833333  0.0833333  0.0833333
        0.238342  0.0799508 0.0446043  0.026387   0.0135135  0.016      0.0272109
        0.282609  0.0799508 0.0446043  0.026387   0.0135135  0.016      0.0272109
    ]
    mb_rate_lb = flatten_group_size(mb_rate_mean[:, :] .* MB_RATE_LB_FACTOR)
    mb_rate_ub = flatten_group_size(mb_rate_mean[:, :] .* MB_RATE_UB_FACTOR)

    coral_params = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Coral)

    lin_ext_idx, mbrate_idx, dhw_tol_mean_idx, dhw_tol_std_idx =
        extract_param_group_idx.([coral_params], target_param_names())

    coral_param_idx = vcat(lin_ext_idx, mbrate_idx, dhw_tol_mean_idx, dhw_tol_std_idx)
    coral_params = coral_params[sort(coral_param_idx), :]

    _, _, dhw_tol_mean_idx, dhw_tol_std_idx =
        extract_param_group_idx.([coral_params], target_param_names())
    lin_ext_idx = sc_fg_param_idxs("linear_extension", coral_params)
    mbrate_idx = sc_fg_param_idxs("mb_rate", coral_params)

    mbrate_scale_idx = findall(occursin.(Ref("mb_rate_scale"), string.(coral_params.fieldname)))
    lin_ext_scale_idx = findall(occursin.(Ref("linear_extension_scale"), string.(coral_params.fieldname)))

    sample_bounds = collect(zip(coral_params.lower_bound, coral_params.upper_bound))
    # Anchor the dist_mean box to the :legacy literature values, not the spec default: the
    # :calib defaults are themselves a calibration output, so anchoring to them would let the
    # search box drift with the ADRIA pin and re-centre on each round's own answer. dist_std
    # has no :legacy variant and needs none - its default is literature-derived.
    dist_mean_anchor = ADRIA.dist_mean(; version=:legacy, n_sizes=N_SIZE_CLASSES)
    sample_bounds[dhw_tol_mean_idx] .= collect(zip(
        dist_mean_anchor .* DHW_TOL_LB_FACTOR,
        dist_mean_anchor .* DHW_TOL_UB_FACTOR
    ))
    sample_bounds[dhw_tol_std_idx] .= collect(zip(
        coral_params[dhw_tol_std_idx, :val] .* DHW_TOL_LB_FACTOR,
        coral_params[dhw_tol_std_idx, :val] .* DHW_TOL_UB_FACTOR
    ))
    set_bounds!(sample_bounds, lin_ext_idx, lin_ext_lb, lin_ext_ub)
    set_bounds!(sample_bounds, mbrate_idx, mb_rate_lb, mb_rate_ub)

    n_mbrate_scale = length(mbrate_scale_idx)
    set_bounds!(sample_bounds, mbrate_scale_idx,
        fill(MB_RATE_SCALE_LB, n_mbrate_scale), fill(MB_RATE_SCALE_UB, n_mbrate_scale))

    n_lin_ext_scale = length(lin_ext_scale_idx)
    set_bounds!(sample_bounds, lin_ext_scale_idx,
        fill(LIN_EXT_SCALE_LB, n_lin_ext_scale), fill(LIN_EXT_SCALE_UB, n_lin_ext_scale))

    coral_start_idx = 1
    coral_end_idx = length(sample_bounds)

    biogroups_ordering = sort(unique(dom.loc_data.CB_CALIB_GROUPS))
    n_biogroups = length(biogroups_ordering)

    growth_acc_start_idx = coral_end_idx + 1
    biogroup_accel_bounds = Matrix{Tuple{Float64,Float64}}(undef, n_biogroups, N_GROWTH_ACCEL_PARAMS)
    biogroup_accel_bounds[:, 1] .= [GROWTH_ACCEL_STEEPNESS_BOUNDS]
    biogroup_accel_bounds[:, 2] .= [GROWTH_ACCEL_HEIGHT_BOUNDS]
    biogroup_accel_bounds[:, 3] .= [GROWTH_ACCEL_MIDPOINT_BOUNDS]
    append!(sample_bounds, accel_params_array_to_vec(biogroup_accel_bounds))
    growth_acc_end_idx = length(sample_bounds)

    sc_dist_start_idx = growth_acc_end_idx + 1
    append!(sample_bounds, fill((SC_DIST_LB, SC_DIST_UB), n_biogroups))
    sc_dist_end_idx = length(sample_bounds)

    param_idxs = [
        coral_start_idx, coral_end_idx,
        growth_acc_start_idx, growth_acc_end_idx,
        sc_dist_start_idx, sc_dist_end_idx
    ]

    coral_param_names = coral_params.fieldname
    growth_accel_names = accel_params_array_to_vec(
        generate_growth_accel_names(collect(1:n_biogroups))
    )

    return CalibConfig(
        sample_bounds, biogroups_ordering, param_idxs, coral_param_names,
        growth_accel_names
    )
end
