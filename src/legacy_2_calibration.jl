# ! This file should probably be deleted

using ADRIA: bleaching_mortality!
using BlackBoxOptim: init_rng!

include("1_setup.jl")

# ---- SETUP SAMPLE BOUNDS -------

# base_location_vector = [
#     (0.0,  0.95), # habitable cover
#     (0.01, 1.0), # taxa weightings
#     (0.01, 1.0),
#     (0.01, 1.0),
#     (0.01, 1.0),
#     (0.01, 1.0),
#     (0.0,  5.0), # size distribution
#     (0.0,  5.0),
#     (0.0,  5.0),
#     (0.0,  5.0),
#     (0.0,  5.0),
# ]

# Define parameter space to scan over
# sample_bounds = repeat(base_location_vector, n_classifications)
sample_bounds::Vector{Tuple{Float64,Float64}} = []
coral_params = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Coral)

taxa_names = [
    "tabular_Acropora",
    "corymbose_Acropora",
    "corymbose_non_Acropora",
    "small_massives",
    "large_massives"
]
szs = 1:7

coral_p_names::Vector{Symbol} = []
coral_bounds::Vector{Tuple{Float64,Float64}} = []

# Create growth bounds
size_widths = ADRIA.bin_widths()
for (t_idx, taxa) in enumerate(taxa_names)
    for s in szs
        push!(
            coral_p_names,
            Symbol(taxa * "_" * string(t_idx) * "_" * string(s) * "_" * "linear_extension")
        )
        push!(coral_bounds, (size_widths[t_idx, s] / 50, size_widths[t_idx, s] * 0.8))
    end
end

# Add mortality bounds
for (t_idx, taxa) in enumerate(taxa_names)
    for s in szs
        push!(
            coral_p_names,
            Symbol(taxa * "_" * string(t_idx) * "_" * string(s) * "_" * "mb_rate")
        )
        if s < 4
            push!(coral_bounds, (0.01, 0.6))
        else
            push!(coral_bounds, (0.01, t_idx > 3 ? 0.15 : 0.4))
        end
    end
end

# Add fecundity bounds
for (t_idx, taxa) in enumerate(taxa_names)
    for s in 4:7
        push!(
            coral_p_names,
            Symbol(taxa * "_" * string(t_idx) * "_" * string(s) * "_" * "fecundity")
        )
        push!(coral_bounds, (Float64(1e4), Float64(1e6)))
    end
end

# Add dist_mean bounds
for (t_idx, taxa) in enumerate(taxa_names)
    for s in 1:7
        tmp_nm = Symbol(taxa * "_" * string(t_idx) * "_" * string(s) * "_" * "dist_mean")
        tmp_val = coral_params[findfirst(x -> x == tmp_nm, coral_params.fieldname), :val]
        push!(coral_p_names, tmp_nm)
        push!(coral_bounds, (tmp_val * 0.9, tmp_val * 1.1))
    end
end

# Add dist_std bounds
for (t_idx, taxa) in enumerate(taxa_names)
    for s in 1:7
        tmp_nm = Symbol(taxa * "_" * string(t_idx) * "_" * string(s) * "_" * "dist_std")
        tmp_val = coral_params[findfirst(x -> x == tmp_nm, coral_params.fieldname), :val]
        push!(coral_p_names, tmp_nm)
        push!(coral_bounds, (tmp_val * 0.9, tmp_val * 1.1))
    end
end

coral_start_idx = length(sample_bounds) + 1
append!(sample_bounds, coral_bounds)
coral_end_idx = length(sample_bounds)

n_taxa = 5
n_limited_locs = length(LIMITED_LOCATIONS)
n_factors = 3 # growth, mortality, fecundity

location_coef = repeat([(0.5, 2.0)], n_taxa * n_limited_locs * n_factors)

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
    err_counts[err_counts.==0] .= 1

    return err_series ./ err_counts
end

"""
Calculate the functional group correlation and temporal correlation between aggregated ltmp
data created for reefmod.
"""
function reef_taxa_error(
    cover;
    rm_ltmp_taxa=rm_ltmp_taxa,
    dom_idxs=dom_idxs
)
    fg_corr::Float64 = 0.0
    for (j, idx) in enumerate(dom_idxs)
        non_missing_mask = (!).(ismissing.(rm_ltmp_taxa[:, 1, j]))
        for id in eachindex(2008:2022)[non_missing_mask]
            fg_corr +=
                cor(cover[id, :, idx], rm_ltmp_taxa[id, :, j]) ./ count(non_missing_mask)
        end
    end
    return fg_corr ./ length(dom_idxs)
end

function insert_init_loc_cover!(
    dom;
    raw_ltmp_reef_data=raw_ltmp_reef_data,
    rm_ltmp_taxa=rm_ltmp_taxa,
    target_dom_idxs=TARGET_DOM_IDXS
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
    reef_error(cover; ltmp_reef_data=ltmp_reef_data)::Vector{Float64}

Calculate the error between ltmp observations and the given cover array.
"""
function reef_error(
    cover;
    raw_ltmp_reef_data=raw_ltmp_reef_data,
    target_dom_idxs=TARGET_DOM_IDXS
)::Vector{Float64}
    err_series::Vector{Float64} = zeros(Float64, 15)
    tmp_err::Vector{Float64} = zeros(Float64, 15)
    err_counts::Vector{Int64} = zeros(Int64, 15)
    not_missing::BitVector = BitVector(repeat([true], 15))
    for (row_idx, ltmp_row) in enumerate(eachrow(raw_ltmp_reef_data))
        if target_dom_idxs[row_idx] == -1
            continue
        end
        not_missing .= (!).(ismissing.(ltmp_row))
        min_arg = argmin(ltmp_row[not_missing])
        max_arg = argmax(ltmp_row[not_missing])
        tmp_err[not_missing] .= MAEE_series(
            cover[not_missing, TARGET_DOM_IDXS[row_idx]], ltmp_row[not_missing]
        )
        tmp_err[not_missing][[min_arg, max_arg]] .*= 2
        err_series[not_missing] .+= tmp_err[not_missing]
        err_counts[not_missing] .+= 1
    end
    if any(err_counts .== 0)
        @debug "No reef level observation data for some years."
        err_counts[err_counts.==0] .= 1
    end

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
        tmp[tmp.<0] .= 0.0
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
        tmp[tmp.<0] .= 0.0
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
    START_YEAR=START_YEAR,
    END_YEAR=END_YEAR,
    coral_param_names=coral_p_names,
    param_idxs=[coral_start_idx, coral_end_idx, loc_coef_start_idx, loc_coef_end_idx],
    loc_idxs=dom_idxs
)
    scen = ADRIA.param_table(dom)
    coral_param_values = init_values[param_idxs[1]:param_idxs[2]]
    scen[1, coral_param_names] = coral_param_values
    corals = ADRIA.to_coral_spec(scen[1, :])

    loc_k_areas = ADRIA.site_k_area(dom)
    loc_areas = ADRIA.loc_area(dom)

    scale_factors::Array{Float64,3} = reshape(
        init_values[param_idxs[3]:param_idxs[4]], (5, 4, 3)
    )
    bleaching_threshold::Float64 = init_values[end]

    linear_ext::Matrix{Float64} = _to_group_size(corals.linear_extension)
    survival_r::Matrix{Float64} = _to_group_size(corals.mb_rate)

    lin_ext_validity = validate_linear_extension_coefficients(
        linear_ext, scale_factors[:, :, 1]
    )
    mortality_validity = validate_mortality_coefficients(
        survival_r, scale_factors[:, :, 2]
    )
    if mortality_validity + lin_ext_validity != 0.0
        return mortality_validity + lin_ext_validity + 6e5
    end

    res = nothing
    try
        res = ADRIA.run_model(dom, scen[1, :], scale_factors, loc_idxs, 4.0)
    catch err
        return 5e5 + sum(linear_ext) + sum(survival_r) + sum(scale_factors)
    end

    # north_mean_cover = zeros(size(res.raw, 1))
    # center_mean_cover = zeros(size(res.raw, 1))
    # south_mean_cover = zeros(size(res.raw, 1))

    # Hacky manual specification of observation years to use to compare against LTMP data
    # which misses some years (ignoring 2023, the last entry)
    # ref_years = START_YEAR:END_YEAR
    # comp_years_north = (ref_years .∈ [ltmp_north[ltmp_north.Year .>= START_YEAR, :Year]])
    # comp_years_center = (ref_years .∈ [ltmp_central[ltmp_central.Year .>= START_YEAR, :Year]])
    # comp_years_south = (ref_years .∈ [ltmp_south[ltmp_south.Year .>= START_YEAR, :Year]])

    # for ts in axes(res.raw, 1)
    # rac = vec(sum(res.raw[ts, :, :], dims=1) .* loc_k_areas') ./ loc_areas
    # north_mean_cover[ts] = mean(rac[NORTH_MASK .&& ltmp_loc_mask])
    # center_mean_cover[ts] = mean(rac[CENTRAL_MASK .&& ltmp_loc_mask])
    # south_mean_cover[ts] = mean(rac[SOUTH_MASK .&& ltmp_loc_mask])
    # end

    # north_mean_cover = north_mean_cover[comp_years_north]
    # center_mean_cover = center_mean_cover[comp_years_center]
    # south_mean_cover = south_mean_cover[comp_years_south]

    # Determine performance
    # north_perf_series = MAEE_series(north_mean_cover, north_obs)
    # central_perf_series = MAEE_series(center_mean_cover, central_obs)
    # south_perf_series = MAEE_series(south_mean_cover, south_obs)

    # north_perf = temporal_variability(north_perf_series)
    # central_perf = temporal_variability(central_perf_series)
    # south_perf = temporal_variability(south_perf_series)

    loc_cover = dropdims(sum(res.raw; dims=2); dims=2) .* loc_k_areas' ./ loc_areas'

    # class_error_series = class_error(loc_cover)
    reef_error_series = reef_error(loc_cover)

    # class_perf = temporal_variability(class_error_series)
    reef_perf = temporal_variability(reef_error_series)

    # Penalize poor performance at lowest trough
    # north_lowest = argmin(north_obs)
    # central_lowest = argmin(central_obs)
    # south_lowest = argmin(south_obs)
    # trough_score = [
    # north_perf_series[north_lowest],
    # central_perf_series[central_lowest],
    # south_perf_series[south_lowest],
    # ]

    # Penalize poor performance at end of time series
    # end_score = [
    # north_perf_series[end],
    # central_perf_series[end],
    # south_perf_series[end],
    # class_error_series[end] * 3,
    # reef_error_series[end] * 10
    # ]

    taxa_cover = dropdims(
        sum(permutedims(reshape(res.raw, (15, 7, 5, 3806)), (1, 3, 2, 4)); dims=3); dims=3
    )
    fg_corr = reef_taxa_error(taxa_cover)

    last_non_zero = findlast(x -> x != 0, reef_error_series)
    first_non_zero = findfirst(x -> x != 0, reef_error_series)

    return 6.0 * (
        reef_perf + reef_error_series[first_non_zero] + reef_error_series[last_non_zero]
    ) + (1 - fg_corr)
end

init_state = deserialize(init_cover_fn)
construct_cover!(dom, init_state, location_classification.consecutive_classification)

insert_init_loc_cover!(dom)

best_score_file = init_guess_fn
if !@isdefined(best_init_state) && isfile(best_score_file)
    best_init_state = deserialize(best_score_file)
    @assert all(first.(sample_bounds) .<= best_init_state .<= last.(sample_bounds)) "Initial state is out of bounds"
else
    best_init_state = nothing
end

function save_results_callback(oc)
    out_fn = "$(OUT_DIR)/coral_p_calib_fixeds.dat"
    best_state = best_candidate(oc)
    serialize(out_fn, best_state)
    @info "Saving"
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
        CallbackInterval=3600
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

out_fn = joinpath(OUT_DIR, "coral_p_calib_fixeddd.dat")

best_fitness(res)
best_init_state = best_candidate(res)

serialize(out_fn, best_init_state)
