"""
Calibrate linear extension, background mortality and fecundity.
Attempt to calibrate location-specific scaling as well.
"""

using ADRIA: bleaching_mortality!
using BlackBoxOptim: init_rng!

include("./1_setup.jl")
include("./param_bounds.jl")

"""
Reshape the growth acceleration parameters from a vector a matrix of shape [n_parameters ⋅ n_locs]
"""
function reshape_growth_accel_parameters(
    params::Vector{Float64};
    n_locs=length(target_dom_idxs)
)::Matrix{Float64}
    return reshape(params, (3, n_locs))
end

function insert_init_loc_cover!(
    dom,
    lambdas;
    raw_ltmp_reef_data=raw_ltmp_reef_data,
    rm_ltmp_taxa=rm_ltmp_taxa,
    target_dom_idxs=target_dom_idxs
)::Nothing
    for (idx, row_idx) in enumerate(target_dom_idxs)
        size_class_props = size_class_distribution(lambdas[idx], ADRIA.bin_edges()[1, :])
        loc_cov =
            rm_ltmp_taxa[2, :, idx] .* size_class_props' ./ sum(rm_ltmp_taxa[2, :, idx])
        tot_cov =
            raw_ltmp_reef_data[
                idx, findfirst(x -> !ismissing(x), raw_ltmp_reef_data[idx, :])
            ] ./ dom.loc_data.k[row_idx]
        dom.init_coral_cover[:, row_idx] .=
            reshape(permutedims(loc_cov, (2, 1)), (35,)) .* tot_cov
    end

    return nothing
end


"""
    validate_linear_extension_coefficients(linear_ext_vals::Matrix{Float64}, linear_ext_coefs::Vector{Float64})::Bool

Check that the sampled linear extension values and location coefficients guarentee values that
do not exceed the size class bin widths.
"""
function validate_linear_extension_coefficients(
    linear_ext_vals::Matrix{Float64}, linear_ext_coefs::Matrix{Float64}
)::Float64
    size_class_bins::Matrix{Float64} = ADRIA.bin_widths()
    n_locs::Int64 = size(linear_ext_coefs, 2)
    dist_from_valid::Float64 = 0.0

    for j in 1:n_locs
        tmp = (linear_ext_vals .* linear_ext_coefs[:, j]) .- size_class_bins
        tmp[tmp .< 0] .= 0.0
        dist_from_valid += sum(tmp)
    end

    return dist_from_valid
end

"""
    validate_mortality_coefficients(mortality_vals::Matrix{Float64}, mortality_coefs::Vector{Float64})::Bool

Ensure that the sampled mortality values and location coefficients guarentee values between
0 and 1.
"""
function validate_mortality_coefficients(
    mortality_vals::Matrix{Float64}, mortality_coefs::Matrix{Float64}
)::Float64
    n_locs::Int64 = size(mortality_coefs, 2)
    survival_vals::Matrix{Float64} = 1 .- mortality_vals

    dist_from_valid::Float64 = 0.0
    for j in 1:n_locs
        tmp = survival_vals .* mortality_coefs[:, j] .- 1
        tmp[tmp .< 0] .= 0.0
        dist_from_valid += sum(tmp)
    end

    return dist_from_valid
end

function _to_group_size(flat_vec::Vector{Float64})::Matrix{Float64}
    return permutedims(reshape(flat_vec, (7, 5)), (2, 1))
end

# For finding minimum index
argmin_missing(x) = argmin(ifelse.(ismissing.(x), Inf, x))

# For finding maximum index
argmax_missing(x) = argmax(ifelse.(ismissing.(x), -Inf, x))

"""
Model simulations end at 2022, so north/central/south obs argument ignores the last entry
which is for 2023.
"""
function obj_func(
    init_values;
    dom_raw=dom,
    location_classification=location_classification.consecutive_classification,
    start_year=start_year,
    end_year=end_year,
    coral_param_names=CORAL_PARAM_NAMES,
    param_idxs=PARAM_IDXS,
    loc_idxs=target_dom_idxs
)
    dom = deepcopy(dom_raw)

    scen = ADRIA.param_table(dom)
    coral_param_values = init_values[param_idxs[1]:param_idxs[2]]
    scen[1, coral_param_names] = coral_param_values

    growth_acc_params = reshape_growth_accel_parameters(
        init_values[param_idxs[5]:param_idxs[6]]
    )

    corals = ADRIA.to_coral_spec(scen[1, :])

    loc_k_areas = ADRIA.site_k_area(dom)
    loc_areas = ADRIA.loc_area(dom)

    # Location-specific scaling for [groups] ⋅ [locations] ⋅ [growth, mortality, fecundity]
    scale_factors::Array{Float64,3} = reshape(
        init_values[param_idxs[3]:param_idxs[4]], (5, 4, 3)
    )

    linear_ext::Matrix{Float64} = _to_group_size(corals.linear_extension)
    survival_r::Matrix{Float64} = _to_group_size(corals.mb_rate)

    lin_ext_validity = validate_linear_extension_coefficients(
        linear_ext, scale_factors[:, :, 1]
    )
    mortality_validity = validate_mortality_coefficients(
        survival_r, scale_factors[:, :, 2]
    )
    if mortality_validity + lin_ext_validity != 0.0
        return 1e6 + (2e6 * (mortality_validity + lin_ext_validity))
    end

    insert_init_loc_cover!(dom, init_values[param_idxs[7]:param_idxs[8]])

    res = nothing
    try
        res = ADRIA.run_model(dom, scen[1, :], scale_factors, growth_acc_params, loc_idxs)
    catch err
        if !(err isa AssertionError)
            rethrow(err)
        elseif contains(sprint(showerror, err), "no method matching")
            # Catch errors unrelated to perturbed parameter values
            rethrow(err)
        end

        # Symmetric error function, where result is 0 if x == y
        # but as y moves further away from x in either positive or negative direction, then
        # the value approaches 1.0 (result is ≈1.0 if the distance between x and y is 10)
        err_func(x, y) = 2.0 / (1.0 + exp(-abs(((x - y) / x) - (y / x)))) - 1.0

        # Determine how far away from the "default" values these are
        comp_params = ADRIA.component_params(dom.model, :Coral)
        coral_fn = comp_params[:, :fieldname]
        lin_ext_idx = contains.(string.(coral_fn), "linear_ext")
        mbrate_idx = contains.(string.(coral_fn), "mb_rate")

        linext_overage =
            err_func.(
                corals.linear_extension, replace(comp_params[lin_ext_idx, :val], 0.0 => 1.0)
            )
        mbrate_overage =
            err_func.(corals.mb_rate, replace(comp_params[mbrate_idx, :val], 0.0 => 1.0))

        return 5e5 + (sum(linext_overage) + sum(mbrate_overage) + sum(scale_factors .> 1.0))
    end

    loc_cover = dropdims(sum(res.raw; dims=2); dims=2) .* loc_k_areas' ./ loc_areas'

    # class_error_series = class_error(loc_cover)
    reef_error_series = reef_error(loc_cover)

    # class_perf = temporal_variability(class_error_series)
    reef_perf = temporal_variability(reef_error_series)

    taxa_cover = dropdims(
        sum(permutedims(reshape(res.raw, (15, 7, 5, 3806)), (1, 3, 2, 4)); dims=3); dims=3
    )
    fg_corr = reef_taxa_error(taxa_cover)

    first_non_zero = findfirst(x -> x != 0, reef_error_series)
    last_non_zero = findlast(x -> x != 0, reef_error_series)

    # Mean Absolute Error for peaks and troughs
    peaks = argmax_missing.(eachrow(raw_ltmp_reef_data))
    troughs = argmin_missing.(eachrow(raw_ltmp_reef_data))
    obs_peaks = raw_ltmp_reef_data[CartesianIndex.([1, 2, 3, 4], peaks)]
    obs_troughs = raw_ltmp_reef_data[CartesianIndex.([1, 2, 3, 4], troughs)]
    modelled_peaks = loc_cover[CartesianIndex.(peaks, target_dom_idxs)]
    modelled_troughs = loc_cover[CartesianIndex.(troughs, target_dom_idxs)]

    peaks_score = mean(abs.(modelled_peaks .- obs_peaks)) * 2.0
    troughs_score = mean(abs.(modelled_troughs .- obs_troughs)) * 2.0

    score = (
        reef_perf +
        reef_error_series[first_non_zero] +
        reef_error_series[last_non_zero] +
        peaks_score +
        troughs_score
    )
    score += fg_corr * 2.0

    return score
end

function restructure_initial_guess!(
    init_guess::Vector{Float64};
    n_locs=length(target_dom_idxs)
)::Vector{Float64}
    # use calibrated growth params as initial guess for location specific
    growth_params::Vector{Float64} = init_guess[(end - 2):end]
    for _ in 1:(n_locs - 1)
        append!(init_guess, growth_params)
    end
    return init_guess
end

init_state = deserialize(init_cover_fn)
construct_cover!(dom, init_state, location_classification.consecutive_classification)

best_score_file = joinpath(OUT_DIR, init_guess_fn)
if isfile(best_score_file)
    best_init_state = deserialize(best_score_file)
    best_init_state = restructure_initial_guess!(best_init_state)
    append!(best_init_state, fill(2.0, length(target_dom_idxs)))
    @assert all(first.(sample_bounds) .<= best_init_state .<= last.(sample_bounds)) "Initial state is out of bounds"
else
    best_init_state = nothing
end

global LAST_SAVE = 0.0

"""
    save_results_callback(
        oc;
        time_interv=1800,
        step_interv=1000,
        result_fn="intermediate_coral_calib.dat"
    )::Nothing

Save intermediate results when either `time_interv` or `step_interv` is reached.
Default values are to save every 30mins or 1000 steps.

Progress plots are created and saved in the format of:

- `calib_progress_[start date and time in UTC]_[number of seconds elapsed].png`
"""
function save_results_callback(
    oc;
    time_interv=1800,
    step_interv=1000,
    result_fn="intermediate_coral_calib.dat"
)::Nothing
    start_time = replace(string(unix2datetime(oc.start_time)), "T"=>"_", ":"=>"")
    elapsed = oc.last_report_time - oc.start_time
    elapsed = round(elapsed; digits=2)

    if LAST_SAVE == 0.0
        # If very first run, save results so we can see how things started
        is_save_point = true
    else
        # Otherwise, check if the required number of steps has elapsed
        is_save_point = oc.num_steps % step_interv == 0

        if !is_save_point
            # If steps have not elapsed, check if the number of seconds has elapsed

            if (oc.last_report_time - LAST_SAVE) > time_interv
                is_save_point = true
            end
        end
    end

    if !is_save_point
        # Save point has not been reached, so exit
        return nothing
    end

    # Otherwise, save intermediate progress!
    global LAST_SAVE = datetime2unix(now(UTC))
    plot_fn = joinpath(OUT_DIR, "calib_progress_$(start_time)_$(elapsed).png")
    calib_fn = joinpath(OUT_DIR, result_fn)
    best_state = best_candidate(oc)
    serialize(calib_fn, best_state)

    interim_res = progress_run(best_state, coral_param_names)
    plot_calibration(interim_res; save_fn=plot_fn)
    @info "Saved intermediate progress"

    return nothing
end

available_threads = Threads.nthreads()

threads_display = available_threads == 1 ? available_threads : available_threads - 1
@info "Using $(threads_display) threads."
if isnothing(best_init_state)
    # Include additional config if using BorgMOEA
    # Method=:borg_moea,
    # FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
    res = bboptimize(
        obj_func;
        SearchRange=sample_bounds,
        MaxSteps=100_000,
        NThreads=available_threads - 1,
        CallbackFunction=save_results_callback,
        CallbackInterval=0  # run at end of every step
    )
else
    @info "Using initial guess."
    res = bboptimize(
        obj_func,
        best_init_state;  # provide an initial solution
        SearchRange=sample_bounds,
        MaxSteps=100_000,
        NThreads=available_threads,
        CallbackFunction=save_results_callback,
        CallbackInterval=0  # run at end of every step
    )
end

out_fn = joinpath(OUT_DIR, "final_calibration_results.dat")

best_fitness(res)
best_init_state = best_candidate(res)

serialize(out_fn, best_init_state)
