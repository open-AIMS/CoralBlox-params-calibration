"""
Calibrate linear extension, background mortality and fecundity.
Attempt to calibrate location-specific scaling as well.
"""

using ADRIA: bleaching_mortality!
using BlackBoxOptim: init_rng!

include("./common/common.jl")
include("./common/cover_construction.jl")

# Define parameter space to scan over
coral_params = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Coral)

# Extract just the target coral parameters
coral_param_idx = extract_param_group_idx(coral_params, "linear_extension")
append!(coral_param_idx, extract_param_group_idx(coral_params, "mb_rate"))
append!(coral_param_idx, extract_param_group_idx(coral_params, "fecundity"))
coral_params = coral_params[sort(coral_param_idx), :]
coral_param_names = coral_params.fieldname

sample_bounds = collect(zip(
    first.(coral_params.dist_params),
    getindex.(coral_params.dist_params, 2)
))

# Adjust bounds for linear extensions
lin_ext_pos = extract_param_group_idx(coral_params, "linear_extension")
size_widths = ADRIA.bin_widths()'[:]  # transpose and flatten

extended_lb = first.(sample_bounds[lin_ext_pos]) .* 0.25
extended_ub = last.(sample_bounds[lin_ext_pos]) .* 2.0
sample_bounds[lin_ext_pos] .= collect(zip(extended_lb, extended_ub))

# Adjust bounds for mortality rate (size: 1 - 3, all groups)
mb_rate_size_1_to_3 = extract_param_group_idx(coral_params, "1_mb_rate")
append!(mb_rate_size_1_to_3, extract_param_group_idx(coral_params, "2_mb_rate"))
append!(mb_rate_size_1_to_3, extract_param_group_idx(coral_params, "3_mb_rate"))
sample_bounds[mb_rate_size_1_to_3] .= fill((0.01, 0.6), length(mb_rate_size_1_to_3))

for grp in 1:5, sz in 1:7
    group_idx = extract_param_group_idx(coral_params, "$(grp)_$(sz)_mb_rate")

    if sz < 4
        sample_bounds[group_idx] .= fill((0.01, 0.6), length(group_idx))
    else
        ub = grp > 3 ? 0.15 : 0.4
        sample_bounds[group_idx] .= fill((0.01, ub), length(group_idx))
    end
end

# Add parameters for location-specific scaling
# Location specific scaling
n_taxa = 5
n_limited_locs = length(limited_locations)
n_factors = 3  # growth, mortality, fecundity

location_coef = fill((0.3, 1.2), n_taxa * n_limited_locs * n_factors)

coral_start_idx = 1
coral_end_idx = length(sample_bounds)

loc_coef_start_idx = coral_end_idx + 1
append!(sample_bounds, location_coef)
loc_coef_end_idx = length(sample_bounds)

"""
    average_class_cover(cover; loc_classes=location_classification.consecutive_classification)::Array{Float64}

Calculate the average cover for each location classification.
"""
function average_class_cover(
    cover;
    loc_classes=location_classification.consecutive_classification
)::Matrix{Float64}
    classes = sort(unique(loc_classes))
    n_tsteps, n_locs = size(cover)
    class_cover::Matrix{Float64} = zeros(Float64, n_tsteps, length(classes))
    class_mask::BitVector = Vector(repeat([true], n_locs))
    for class in classes
        class_mask .= loc_classes .== class
        class_cover[:, class] .= dropdims(mean(cover[:, class_mask]; dims=2); dims=2)
    end

    return class_cover
end

function insert_init_loc_cover!(
    dom;
    raw_ltmp_reef_data=raw_ltmp_reef_data,
    rm_ltmp_taxa=rm_ltmp_taxa,
    target_dom_idxs=target_dom_idxs
)::Nothing
    size_class_props = size_class_distribution(2.0, ADRIA.bin_edges()[1, :])
    for (idx, row_idx) in enumerate(target_dom_idxs)
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
Calculate the functional group correlation and temporal correlation between aggregated ltmp
data created for reefmod.

Uses the complement of the absolute pearson correlation coefficient such that 0 indicates
a perfect fit, and 1 indicates no correlation.
"""
function reef_taxa_error(
    cover;
    rm_ltmp_taxa=rm_ltmp_taxa,
    dom_idxs=target_dom_idxs
)
    fg_corr::Float64 = 0.0
    for (j, idx) in enumerate(dom_idxs)
        non_missing_mask = (!).(ismissing.(rm_ltmp_taxa[:, 1, j]))
        for id in eachindex(2008:2022)[non_missing_mask]
            fg_corr +=
                abs(cor(cover[id, :, idx], rm_ltmp_taxa[id, :, j])) ./ count(non_missing_mask)
        end
    end

    return 1.0 - (fg_corr ./ length(dom_idxs))
end

"""
    reef_error(cover; ltmp_reef_data=ltmp_reef_data)::Vector{Float64}

Calculate the error between ltmp observations and the given cover array.
"""
function reef_error(
    cover;
    ltmp_obs=raw_ltmp_reef_data,
    target_dom_idxs=target_dom_idxs
)::Vector{Float64}
    err_series::Vector{Float64} = zeros(Float64, 15)
    tmp_err::Vector{Float64} = zeros(Float64, 15)
    err_counts::Vector{Int64} = zeros(Int64, 15)
    not_missing::BitVector = BitVector(fill(true, 15))
    for (row_idx, loc_obs) in enumerate(eachrow(ltmp_obs))
        if target_dom_idxs[row_idx] == -1
            continue
        end

        not_missing .= (!).(ismissing.(loc_obs))
        min_arg = argmin(loc_obs[not_missing])
        max_arg = argmax(loc_obs[not_missing])
        tmp_err[not_missing] .= MAEE_series(
            cover[not_missing, target_dom_idxs[row_idx]], loc_obs[not_missing]
        )

        # Apply double the weight on the peak/trough of the time series
        tmp_err[not_missing][[min_arg, max_arg]] .*= 2.0
        err_series[not_missing] .+= tmp_err[not_missing]
        err_counts[not_missing] .+= 1.0
    end

    if any(err_counts .== 0)
        @debug "No reef level observation data for some years."
        err_counts[err_counts .== 0] .= 1
    end

    return err_series ./ err_counts
end

"""
    class_error(cover; obs_class_data=manta_tow_classes)::Vector{Float64}

Calculate the average class level error per year.
"""
function class_error(
    cover;
    manta_tow_mean=manta_tow_mean,
    manta_tow_std=manta_tow_std
)::Vector{Float64}
    # Preallocations
    err_series::Vector{Float64} = zeros(Float64, 15)
    err_counts::Vector{Int64} = zeros(Int64, 15)
    not_missing::BitVector = BitVector(repeat([true], 15))

    # Dims ~ [timesteps ⋅ classes]
    class_cover::Matrix{Float64} = average_class_cover(cover)

    for (idx, class) in enumerate(manta_tow_mean.class)
        if class == -1
            continue
        end
        not_missing .= (!).(ismissing.(manta_tow_mean[idx, :]))
        err_series[not_missing] .+=
            abs.(
                (
                    manta_tow_mean[idx, not_missing] .- class_cover[not_missing, class]
                ) ./ manta_tow_std[idx, not_missing]
            )
        err_counts[not_missing] .+= 1
    end
    err_counts[err_counts .== 0] .= 1

    return err_series ./ err_counts
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

"""
Model simulations end at 2022, so north/central/south obs argument ignores the last entry
which is for 2023.
"""
function obj_func(
    init_values;
    dom=dom,
    location_classification=location_classification.consecutive_classification,
    start_year=start_year,
    end_year=end_year,
    coral_param_names=coral_param_names,
    param_idxs=[coral_start_idx, coral_end_idx, loc_coef_start_idx, loc_coef_end_idx],
    loc_idxs=target_dom_idxs
)
    scen = ADRIA.param_table(dom)
    coral_param_values = init_values[param_idxs[1]:param_idxs[2]]
    scen[1, coral_param_names] = coral_param_values
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

    res = nothing
    try
        res = ADRIA.run_model(dom, scen[1, :], scale_factors, loc_idxs)
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

        linext_overage = err_func.(corals.linear_extension, replace(comp_params[lin_ext_idx, :val], 0.0 => 1.0))
        mbrate_overage = err_func.(corals.mb_rate, replace(comp_params[mbrate_idx, :val], 0.0 => 1.0))

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

    score = reef_perf + reef_error_series[first_non_zero] + reef_error_series[last_non_zero]
    score += fg_corr * 2.0

    return score
end

init_state = deserialize(init_cover_fn)
construct_cover!(dom, init_state, location_classification.consecutive_classification)

insert_init_loc_cover!(dom)

best_score_file = joinpath(OUT_DIR, init_guess_fn)
if !@isdefined(best_init_state) && isfile(best_score_file)
    best_init_state = deserialize(best_score_file)
    @assert all(first.(sample_bounds) .<= best_init_state .<= last.(sample_bounds)) "Initial state is out of bounds"
else
    best_init_state = nothing
end

"""
    save_results_callback(
        oc;
        time_interv=3600,
        step_interv=10_000,
        result_fn="intermediate_coral_calib.dat"
    )::Nothing

Save intermediate results when at least one of the time or step interval is hit.
"""
function save_results_callback(
    oc;
    time_interv=3600,
    step_interv=10_000,
    result_fn="intermediate_coral_calib.dat"
)::Nothing
    is_save_point = oc.num_steps % step_interv == 0
    elapsed = (oc.last_report_time - oc.start_time)
    is_save_time = elapsed > 0 ? elapsed % time_interv == 0 : false

    if !is_save_point && !is_save_time
        return nothing
    end

    out_fn = joinpath(OUT_DIR, result_fn)
    best_state = best_candidate(oc)
    serialize(out_fn, best_state)
    @info "Saved intermediate"
    return nothing
end

@info "Using $(Threads.nthreads()-1) threads."
if !@isdefined(best_init_state) || isnothing(best_init_state)
    # Include additional config if using BorgMOEA
    # Method=:borg_moea,
    # FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
    res = bboptimize(
        obj_func;
        SearchRange=sample_bounds,
        MaxSteps=100_000,
        NThreads=Threads.nthreads() - 1,
        CallbackFunction=save_results_callback,
        CallbackInterval=0  # run at end of every step
    )
elseif !isnothing(best_init_state)
    @info "Using initial guess."
    res = bboptimize(
        obj_func,
        best_init_state;  # provide an initial solution
        SearchRange=sample_bounds,
        MaxSteps=100_000,
        NThreads=Threads.nthreads() - 1
    )
end

out_fn = joinpath(OUT_DIR, "coral_p_calib_last.dat")

best_fitness(res)
best_init_state = best_candidate(res)

serialize(out_fn, best_init_state)
