"""
Methods to run and plot results.
"""

using Dates

"""
    calib_run(dom, params, coral_param_names, param_idxs, loc_idxs)

Run ADRIA-CoralBlox with the provided parameters.

# Arguments
- `dom` : ADRIA Domain
- `params` : Calibrated parameter values
- `coral_param_names` : Parameter names
- `param_idxs` : Indices of parameter sets in parameter set
                 (start of coral, end of coral, start of location scaling, end of location scaling)
- `loc_idxs` : Indices of target locations in reef list

# Returns
ADRIA run results
"""
function calib_run(dom_raw, params, coral_param_names, param_idxs, loc_idxs)
    dom = deepcopy(dom_raw)

    scale_factors::Array{Float64,3} = reshape(
        params[param_idxs[3]:param_idxs[4]], (5, 4, 3)
    )

    scen = ADRIA.param_table(dom)

    growth_acc_params = reshape_growth_accel_parameters(params[param_idxs[5]:param_idxs[6]])

    insert_init_loc_cover!(dom, params[param_idxs[7]:param_idxs[8]])

    coral_param_values = params[param_idxs[1]:param_idxs[2]]
    scen[1, coral_param_names] = coral_param_values

    return ADRIA.run_model(dom, scen[1, :], scale_factors, growth_acc_params, loc_idxs)
end

"""
    convert_to_ltmp_values(res)

Transform ADRIA-CoralBlox results to allow comparison with LTMP data.
"""
function convert_to_ltmp_values(res)
    return dropdims(sum(res.raw, dims=2), dims=2)
end

"""
    progress_run(interm_params, coral_param_names)

Run ADRIA-CoralBlox with calibrated parameters.
"""
function progress_run(interm_params, coral_param_names)
    coral_start_idx = 1
    coral_end_idx = length(coral_param_names)

    loc_coef_start_idx = coral_end_idx + 1
    loc_coef_end_idx = loc_coef_start_idx + 5 * 4 * 3 - 1

    growth_start_idx = loc_coef_end_idx + 1
    growth_end_idx = growth_start_idx + 3 * 4 - 1

    sc_dist_start_idx = growth_end_idx + 1
    sc_dist_end_idx = length(interm_params)

    calib_res = calib_run(
        dom,
        interm_params,
        coral_param_names,
        [
            coral_start_idx, coral_end_idx,
            loc_coef_start_idx, loc_coef_end_idx,
            growth_start_idx, growth_end_idx,
            sc_dist_start_idx, sc_dist_end_idx
        ],
        target_dom_idxs
    )

    return calib_res
end

"""
    plot_calibration(calib_res; save_fn="calib_progress.png")

Create a plot of the four locations targeted for calibration.
"""
function plot_calibration(calib_res; save_fn="calib_progress.png")
    modelled_locs = convert_to_ltmp_values(calib_res)[:, target_dom_idxs]
    obs_locs = raw_ltmp_reef_data'

    f = Figure(size=(900, 600))
    max_col = 2
    row = Ref(1)
    col = Ref(1)
    for (i, loc) in enumerate(target_dom_idxs)
        reef_name = dom.loc_data[loc, :GBR_NAME]
        reef_id = dom.loc_data[loc, :UNIQUE_ID]

        ltmp_loc_pos = findfirst(x->!ismissing(x) && x==reef_id, ltmp_reef_data.RME_UNIQUE_ID)

        rmse_, benchmark_, cc_, maee_, bias_ = collect_error_stats(calib_res.raw, ltmp_loc_pos)

        rmse_ = trunc(rmse_, digits=4)
        benchmark_ = trunc(benchmark_, digits=4)
        cc_ = trunc(cc_, digits=4)
        maee_ = trunc(maee_, digits=4)
        bias_ = trunc(bias_, digits=4)

        err_report_str = "RMSE: $(rmse_) | μ bnch: $(benchmark_) | PCC: $(cc_) | MAEE: $(maee_) | BIAS: $(bias_)"
        title_text = rich("$reef_name\n$(reef_id)\n", rich(err_report_str, fontsize=9))
        ax = Axis(
            f[row[], col[]],
            title=title_text,
            xticks=(1:15, string.(start_year:end_year)),
            xticklabelrotation=45
        )
        scatter!(ax, obs_locs[:, i], color=(:red, 0.5))
        lines!(ax, modelled_locs[:, i])

        if col[] < max_col
            col[] += 1
        else
            row[] += 1
            col[] = 1
        end
    end

    linkaxes!(filter(x -> x isa Axis, f.content)...)

    save(save_fn, f)

    return nothing
end
