"""
Methods to run and plot results.
"""

using Dates

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
function progress_run(interm_params)

    new_dom, new_scen, growth_acc_params, scale_factors = setup_run(dom, interm_params)

    calib_res = ADRIA.run_model(
        new_dom,
        new_scen[1, :],
        scale_factors,
        growth_acc_params,
        target_dom_idxs
    )

    return calib_res
end

"""
    plot_calibration(calib_res; save_fn="calib_progress.png")

Create a plot of the four locations targeted for calibration.
"""
function plot_calibration(calib_res; save_fn="calib_progress.png")

    f = Figure(size=(900, 600))
    max_col = 2
    row = Ref(1)
    col = Ref(1)
    for (i, loc) in enumerate(target_dom_idxs)
        reef_id = dom.loc_data[loc, :UNIQUE_ID]

        ltmp_loc_pos = findfirst(x->!ismissing(x) && x==reef_id, ltmp_reef_data.RME_UNIQUE_ID)

        title_text = construct_location_err_title(calib_res.raw, ltmp_loc_pos)
        Axis(
            f[row[], col[]],
            title=title_text,
            xticks=(1:15, string.(start_year:end_year)),
            xticklabelrotation=45
        )
        plot_modelled_v_ltmp(calib_res.raw, ltmp_loc_pos, loc)

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
