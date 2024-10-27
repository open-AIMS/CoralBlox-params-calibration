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
function calib_run(dom, params, coral_param_names, param_idxs, loc_idxs)
    scale_factors::Array{Float64,3} = reshape(
        params[param_idxs[3]:param_idxs[4]], (5, 4, 3)
    )

    scen = ADRIA.param_table(dom)
    coral_param_values = params[param_idxs[1]:param_idxs[2]]
    scen[1, coral_param_names] = coral_param_values

    return ADRIA.run_model(dom, scen[1, :], scale_factors, loc_idxs)
end

"""
    convert_to_ltmp_values(res)

Transform ADRIA-CoralBlox results to allow comparison with LTMP data.
"""
function convert_to_ltmp_values(res)
    return dropdims(sum(res.raw, dims=2), dims=2)
end

function plot_calibration(param_filepath, coral_param_names)
    interm = deserialize(param_filepath)

    coral_start_idx = 1
    coral_end_idx = length(coral_param_names)

    loc_coef_start_idx = coral_end_idx + 1
    loc_coef_end_idx = length(interm)

    calib_res = calib_run(
        dom,
        interm,
        coral_param_names,
        [coral_start_idx, coral_end_idx, loc_coef_start_idx, loc_coef_end_idx],
        target_dom_idxs
    )

    modelled_locs = convert_to_ltmp_values(calib_res)[:, target_dom_idxs]
    obs_locs = raw_ltmp_reef_data'

    f = Figure(size=(900, 600))
    max_col = 2
    row = Ref(1)
    col = Ref(1)
    for (i, loc) in enumerate(target_dom_idxs)
        reef_name = dom.loc_data[loc, :GBR_NAME]
        reef_id = dom.loc_data[loc, :UNIQUE_ID]

        ax = Axis(
            f[row[], col[]],
            title="$reef_name\n$(reef_id)",
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

    right_now = replace(string(now()), ":"=>"", "."=>"_")
    save(joinpath(OUT_DIR, "calib_progress_$(right_now).png"), f)

    return nothing
end

# plot_calibration(joinpath(OUT_DIR, "coral_p_calib_last.dat"), coral_param_names)
