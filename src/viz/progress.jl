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
    progress_run(
        interm_params, dom, param_names, growth_accel_names, depth_atten_names,
        dist_std_names, dist_mean_names, param_idxs, observations, biogroup_ord
    )

Run ADRIA-CoralBlox with calibrated parameters.
"""
function progress_run(
    interm_params,
    dom,
    param_names::Vector{Symbol},
    growth_accel_names::Vector{String},
    depth_atten_names::Vector{String},
    dist_std_names::Vector{String},
    dist_mean_names::Vector{String},
    param_idxs::Vector{Int64},
    observations::LocationDataStore,
    biogroup_ord::Vector{Int64}
)
    new_dom, new_scen = setup_run(
        dom,
        interm_params;
        param_names=param_names,
        growth_accel_names=growth_accel_names,
        depth_atten_names=depth_atten_names,
        dist_std_names=dist_std_names,
        dist_mean_names=dist_mean_names,
        param_idxs=param_idxs,
        observations=observations,
        biogroup_ord=biogroup_ord,
    )

    calib_res = ADRIA.run_model(
        new_dom,
        new_scen[1, :];
        apply_allee_effect=false,
    )

    return calib_res
end
