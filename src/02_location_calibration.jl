"""
Calibrate linear extension, background mortality and fecundity.
Attempt to calibrate location-specific scaling as well.
"""

using ADRIA: bleaching_mortality!
using BlackBoxOptim

include("plot/plot.jl")
include("./common/param_bounds.jl")

"""
    validate_gbr_wide_scalar_mean(linear_ext_scalar::Matrix{Float64}, survival_scalar::Matrix{Float64})::Float64

The average parameter scalar should between [0.95, 1.05] for linear extension and
[-0.05, 0.05] for background mortality.
"""
function validate_gbr_wide_scalar_mean(
    linear_ext_scalar::Matrix{Float64},
    survival_scalar::Matrix{Float64}
)::Float64
    mean_lin_ext_scalar::Vector{Float64} = dropdims(mean(linear_ext_scalar, dims=2), dims=2)
    mean_survival_scalar::Vector{Float64} = dropdims(mean(survival_scalar, dims=2), dims=2)

    lin_ext_within_bounds::BitVector = (!).(0.95 .<= mean_lin_ext_scalar .<= 1.05)
    survival_within_bounds::BitVector = (!).(-0.05 .<= mean_survival_scalar .<= 0.05)
    return sum(abs.(mean_lin_ext_scalar .- 1.0)[lin_ext_within_bounds]) +
           sum(abs.(mean_survival_scalar)[survival_within_bounds])
end


"""
    validate_linear_extension_coefficients(linear_ext_vals::Matrix{Float64}, linear_ext_coefs::Vector{Float64})::Bool

Check that the sampled linear extension values and location coefficients guarantee values
that do not exceed the size class bin widths.
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

function _to_group_size(flat_vec::Vector{Float64})::Matrix{Float64}
    return permutedims(reshape(flat_vec, (N_SIZE_CLASSES, N_TAXA)), (2, 1))
end

# For finding minimum index
argmin_missing(x) = argmin(ifelse.(ismissing.(x), Inf, x))

# For finding maximum index
argmax_missing(x) = argmax(ifelse.(ismissing.(x), -Inf, x))

# Relative deviation of `calibrated` from `default`, mapped to [0, 1) via a sigmoid.
# Zero when calibrated == default; approaches 1 as |calibrated - default| / |default| grows.
# When default is 0 the denominator is set to 1 to avoid divide-by-zero.
function _relative_deviation_penalty(calibrated::Float64, default::Float64)::Float64
    denom = default == 0.0 ? 1.0 : abs(default)
    rel_dev = abs(calibrated - default) / denom
    return 2.0 / (1.0 + exp(-rel_dev)) - 1.0
end

"""
Model simulations end at 2022, so north/central/south obs argument ignores the last entry
which is for 2023.
"""
function obj_func(
    init_values;
    dom_raw=dom,
    location_classification=location_classification.consecutive_classification,
    coral_param_names=CORAL_PARAM_NAMES,
    param_idxs=PARAM_IDXS,
    observations::LocationDataStore=CALIBRATION_STORE
)
    dom, scen = setup_run(
        dom_raw,
        init_values;
        param_names=coral_param_names,
        param_idxs=param_idxs,
    )

    scale_factors = get_scale_factors(dom, scen)

    corals = ADRIA.to_coral_spec(scen[1, :])

    loc_k_areas = ADRIA.site_k_area(dom)
    loc_areas = ADRIA.loc_area(dom)

    linear_ext::Matrix{Float64} = _to_group_size(corals.linear_extension)

    # If either the linear_extension or mortality are not valid return a proportionally
    # big error value
    gbr_wide_scalar_validity = validate_gbr_wide_scalar_mean(
        scale_factors[:, 1, :], scale_factors[:, 2, :]
    )
    lin_ext_validity = validate_linear_extension_coefficients(
        linear_ext, scale_factors[:, 1, :]
    )
    validity_violation = lin_ext_validity + gbr_wide_scalar_validity
    if validity_violation > 0.0
        # Cannot run the model with invalid parameters; return early with a penalty that
        # scales linearly with the violation magnitude. A flat constant (ignoring magnitude)
        # would make all infeasible candidates look equally bad to the optimizer, blocking
        # directional search toward the feasibility boundary. The floor of 1e4 sits well
        # above any plausible valid score (~25 worst-case) so infeasible candidates always
        # lose to feasible ones in the minimizer.
        return 1e4 + 1e4 * validity_violation
    end

    res = nothing
    try
        res = ADRIA.run_model(dom, scen[1, :])
    catch err
        # Only swallow AssertionErrors known to be triggered by biologically impossible
        # parameter values: (1) recruits_scale_factor <= 0 when linear extension pushes
        # cover past the threshold cap (growth.jl), and (2) survival rate > 1 from extreme
        # mortality parameters (scenario.jl). All other AssertionErrors must propagate —
        # silently converting infrastructure failures to a finite score (~5e5) lets the
        # optimizer treat crash-inducing regions as merely "moderately bad" candidates.
        _param_assert_patterns = (
            "!any(recruits_scale_factor .<= 0)",
            "Survival rate should be <= 1",
        )
        err_msg = sprint(showerror, err)
        is_param_assert = err isa AssertionError &&
            any(contains(err_msg, p) for p in _param_assert_patterns)
        if !is_param_assert
            rethrow(err)
        end

        comp_params = ADRIA.component_params(dom.model, :Coral)
        coral_fn = comp_params[:, :fieldname]
        lin_ext_idx = contains.(string.(coral_fn), "linear_ext")
        mbrate_idx = contains.(string.(coral_fn), "mb_rate")

        lin_ext_overage = _relative_deviation_penalty.(
            corals.linear_extension, comp_params[lin_ext_idx, :val]
        )
        mb_rate_overage = _relative_deviation_penalty.(
            corals.mb_rate, comp_params[mbrate_idx, :val]
        )

        return 5e5 + (sum(lin_ext_overage) + sum(mb_rate_overage) + sum(scale_factors .> 1.0))
    end

    loc_cover = dropdims(sum(res.raw; dims=(2, 3)); dims=(2, 3)) .* loc_k_areas' ./ loc_areas'

    # class_error_series = class_error(loc_cover)
    reef_error_series = reef_error(loc_cover)

    # class_perf = temporal_variability(class_error_series)
    reef_perf = temporal_variability(reef_error_series)

    taxa_cover = dropdims(
        sum(res.raw; dims=3); dims=3
    )
    fg_corr = reef_taxa_error(taxa_cover)

    first_non_zero = findfirst(x -> x != 0, reef_error_series)
    last_non_zero = findlast(x -> x != 0, reef_error_series)

    # Restrict to survey locations with at least one valid observation. An all-missing row
    # causes argmax_missing to return an arbitrary index whose stored value is still Missing,
    # which would silently propagate Missing through to the score.
    has_obs = vec(any(.!ismissing.(observations.ltmp_coral_cover), dims=2))
    valid_survey_rows = findall(has_obs)
    valid_domain_rows = observations.ltmp_cover_to_domain[valid_survey_rows]

    # argmax/argmin_missing return column (year) indices into the cover array, which equal
    # row indices into loc_cover — both share the same [1, n_years] time axis.
    each_obs_reef = eachrow(observations.ltmp_coral_cover[valid_survey_rows, :])
    obs_peak_col = argmax_missing.(each_obs_reef)
    obs_trough_col = argmin_missing.(each_obs_reef)
    obs_peaks = Float64.(observations.ltmp_coral_cover[CartesianIndex.(valid_survey_rows, obs_peak_col)])
    obs_troughs = Float64.(observations.ltmp_coral_cover[CartesianIndex.(valid_survey_rows, obs_trough_col)])
    modelled_peaks = loc_cover[CartesianIndex.(obs_peak_col, valid_domain_rows)]
    modelled_troughs = loc_cover[CartesianIndex.(obs_trough_col, valid_domain_rows)]

    peaks_score = mean(abs.(modelled_peaks .- obs_peaks)) * 0.5
    troughs_score = mean(abs.(modelled_troughs .- obs_troughs)) * 0.5

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

init_state = deserialize(INIT_COVER_PATH)
construct_cover!(dom, init_state, location_classification.consecutive_classification)

best_score_file = joinpath(OUT_DIR, INIT_GUESS_PATH)
if isfile(best_score_file)
    best_init_state = deserialize(best_score_file)
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
    start_time = replace(string(unix2datetime(oc.start_time)), "T" => "_", ":" => "")
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
    region_plot_fn = joinpath(OUT_DIR, "region_plots", "calib_progress_region_$(start_time)_$(elapsed).png")
    taxa_cover_plot_fn = joinpath(OUT_DIR, "taxa_cover", "calib_progress_taxa_cover_$(start_time)_$(elapsed).png")
    taxa_pop_plot_fn = joinpath(OUT_DIR, "taxa_pop", "calib_progress_pop_cover_$(start_time)_$(elapsed).png")
    calib_fn = joinpath(OUT_DIR, result_fn)
    best_state = best_candidate(oc)
    serialize(calib_fn, best_state)

    interim_res = progress_run(best_state)
    f = plot_all_regions(dom, interim_res)
    save(region_plot_fn, f)
    f = taxa_cover_proportions(interim_res.raw)
    save(taxa_cover_plot_fn, f)
    f = taxa_population_proportions(interim_res.raw)
    save(taxa_pop_plot_fn, f)
    @info "Saved intermediate progress"

    return nothing
end

available_threads = Threads.nthreads()

mkpath(joinpath(OUT_DIR, "region_plots"))
mkpath(joinpath(OUT_DIR, "taxa_cover"))
mkpath(joinpath(OUT_DIR, "taxa_pop"))

threads_display = available_threads == 1 ? available_threads : available_threads - 1
@info "Using $(threads_display) threads."

@info "Using RNG seed $(RNG_SEED) for this calibration run (config.toml: operation.rng_seed)."
# NOTE: with NThreads > 0, evaluation/archive-update order depends on OS thread
# scheduling, not just the RNG. Reproducibility here is "deterministic modulo
# thread scheduling" (same seed + hardware + thread count is likely, but not
# guaranteed, bit-exact). Only NThreads=0 (single-threaded) would be bit-exact.
if isnothing(best_init_state)
    # Include additional config if using BorgMOEA
    # Method=:borg_moea,
    # FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
    res = bboptimize(
        obj_func;
        SearchRange=sample_bounds,
        MaxSteps=1_000_000,
        NThreads=available_threads - 1,
        CallbackFunction=save_results_callback,
        CallbackInterval=0,  # run at end of every step
        RngSeed=RNG_SEED,
        RandomizeRngSeed=false
    )
else
    @info "Using initial guess."
    res = bboptimize(
        obj_func,
        best_init_state;  # provide an initial solution
        SearchRange=sample_bounds,
        MaxSteps=1_000_000,
        NThreads=available_threads - 1,
        CallbackFunction=save_results_callback,
        CallbackInterval=0,  # run at end of every step
        RngSeed=RNG_SEED,
        RandomizeRngSeed=false
    )
end

out_fn = joinpath(OUT_DIR, RESULT_FN)

best_fitness(res)
best_init_state = best_candidate(res)

serialize(out_fn, best_init_state)
