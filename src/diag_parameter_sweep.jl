"""
Script to perform a parameter sweep.
"""

using ADRIA: bleaching_mortality!
using BlackBoxOptim: init_rng!

include("1_setup.jl")

target_param_group = "mean_colony_diameter_m"  # change this to target coral factor groups

fig_dir = joinpath(OUT_DIR, "diag_sweep")
try
    mkdir(fig_dir)
catch err
    if !(err isa Base.IOError)
        rethrow(err)
    end
end

function sweep_run(dom, params, coral_param_names, loc_idxs)
    # No location-based scaling
    scale_factors = ones(5, 4, 3)

    scen = ADRIA.param_table(dom)
    scen[1, coral_param_names] .= params

    return ADRIA.run_model(dom, scen[1, :], scale_factors, loc_idxs)
end

# Define parameter space to scan over
coral_params = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Coral)

# Extract just the target coral parameters
target_param_pos = extract_param_group_idx(coral_params, target_param_group)
coral_params = coral_params[target_param_pos, :]
coral_param_names = coral_params.fieldname

# Extract factor range to sweep over
sample_bounds = collect(zip(
    first.(coral_params.dist_params),
    getindex.(coral_params.dist_params, 2) * 1.7  # Note the larger bound here!
))

# Calculate sweeps as five equi-distant steps
step_size = 4   # minimum of four (resolves to five), otherwise there are range issues
steps = (last.(sample_bounds) .- first.(sample_bounds)) ./ step_size
sweep = [ed != 0 ? collect(st:step:ed) : zeros(step_size+1) for (st, step, ed) in zip(first.(sample_bounds), steps, last.(sample_bounds))]

# Collate all parameter values for sweep into a matrix
# [sweep ⋅ factors]
sweep_values = hcat(sweep...)
for (sw_id, pset) in enumerate(eachrow(sweep_values))
    @info "Running sweep $(sw_id) - $(target_param_group)"
    sweep_res = sweep_run(dom, pset, coral_param_names, target_dom_idxs)

    modelled_locs = convert_to_ltmp_values(sweep_res)[:, target_dom_idxs]
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

    save(joinpath(fig_dir, "sweep_$(target_param_group)_$(sw_id).png"), f)
end