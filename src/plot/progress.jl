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

    new_dom, new_scen = setup_run(dom, interm_params)

    calib_res = ADRIA.run_model(
        new_dom,
        new_scen[1, :],
    )

    return calib_res
end

