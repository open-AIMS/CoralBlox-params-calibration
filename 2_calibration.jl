using BlackBoxOptim: init_rng!
include("common.jl")
include("cover_construction.jl")


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
        class_cover[:, class] .= dropdims(mean(cover[:, class_mask], dims=2), dims=2)
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
        err_series[not_missing] .+= abs.((
            manta_tow_mean[idx, not_missing] .- class_cover[not_missing, class]
        ) ./ manta_tow_std[idx, not_missing])
        err_counts[not_missing] .+= 1
    end
    err_counts[err_counts .== 0] .= 1

    return err_series ./ err_counts
end

"""
    reef_error(cover; ltmp_reef_data=ltmp_reef_data)::Vector{Float64}

Calculate the error between ltmp observations and the given cover array.
"""
function reef_error(
    cover;
    raw_ltmp_reef_data=raw_ltmp_reef_data,
    ltmp_reefmod_idxs=ltmp_reefmod_idxs
)::Vector{Float64}
    err_series::Vector{Float64} = zeros(Float64, 15)
    err_counts::Vector{Int64} = zeros(Int64, 15)
    not_missing::BitVector = BitVector(repeat([true], 15))
    for (row_idx, ltmp_row) in enumerate(eachrow(raw_ltmp_reef_data))
        if ltmp_reefmod_idxs[row_idx] == -1
            continue
        end
        not_missing .= (!).(ismissing.(ltmp_row))
        err_series[not_missing] .+= MAEE_series(
            cover[not_missing, ltmp_reefmod_idxs[row_idx]], ltmp_row[not_missing]
        )
        err_counts[not_missing] .+= 1
    end
    if any(err_counts .== 0)
        @warn "No reef level observation data for some years."
        err_counts[err_counts .== 0] .= 1
    end

    return err_series ./ err_counts
end

"""
Model simulations end at 2022, so north/central/south obs argument ignores the last entry
which is for 2023.
"""
function obj_func(
    init_values;
    dom=dom,
    north_cover=ltmp_north[ltmp_north_period, [:lower, :response, :upper]],
    central_cover=ltmp_central[ltmp_central_period, [:lower, :response, :upper]],
    south_cover=ltmp_south[ltmp_south_period, [:lower, :response, :upper]],
    location_classification=location_classification.consecutive_classification,
    n_loc_clusteres=n_classifications,
    start_year=start_year,
    end_year=end_year,
)
    n_cover_start = 11 * n_loc_clusteres
    construct_cover!(dom, init_values, location_classification)
    gen_init_cover = dom.init_coral_cover.data

    loc_k_areas = site_k_area(dom)
    loc_areas = site_area(dom)

    # Check constraints
    total_cover_gen = vec(sum(gen_init_cover, dims=1))

    north_obs = north_cover.response
    central_obs = central_cover.response
    south_obs = south_cover.response

    # If start condition is out of bounds
    start_score = zeros(3)
    total_north = (total_cover_gen[NORTH_MASK] .* loc_k_areas[NORTH_MASK]) ./ loc_areas[NORTH_MASK]
    start_score[1] = abs(mean(total_north) - north_obs[1])
    if !(north_cover.lower[1] .<= mean(total_north) .<= north_cover.upper[1])
        start_score[1] *= 100
        start_score[1] += 1e6
    end

    total_central = (total_cover_gen[CENTRAL_MASK] .* loc_k_areas[CENTRAL_MASK]) ./ loc_areas[CENTRAL_MASK]
    start_score[2] = abs(mean(total_central) - central_obs[1])
    if !(central_cover.lower[1] .<= mean(total_central) .<= central_cover.upper[1])
        start_score[2] *= 100
        start_score[2] += 1e6
    end

    total_south = (total_cover_gen[SOUTH_MASK] .* loc_k_areas[SOUTH_MASK]) ./ loc_areas[SOUTH_MASK]
    start_score[3] = abs(mean(total_south) - south_obs[1])
    if !(south_cover.lower[1] .<= mean(total_south) .<= south_cover.upper[1])
        start_score[3] *= 100
        start_score[3] += 1e6
    end

    if sum(start_score) >= 1e6
        return sum(start_score)
    end

    # If initial cover pass constraint rules, assign them and run
    dom2 = deepcopy(dom)
    dom2.init_coral_cover.data .= gen_init_cover

    scen = ADRIA.param_table(dom2)

    res = nothing
    try
        res = ADRIA.run_model(dom2, scen[1, :])
    catch err
        return sum(start_score) + 5e5
    end

    north_mean_cover = zeros(size(res.raw, 1))
    center_mean_cover = zeros(size(res.raw, 1))
    south_mean_cover = zeros(size(res.raw, 1))

    # Hacky manual specification of observation years to use to compare against LTMP data
    # which misses some years (ignoring 2023, the last entry)
    ref_years = start_year:end_year
    comp_years_north = (ref_years .∈ [ltmp_north[ltmp_north.Year .>= start_year, :Year]])
    comp_years_center = (ref_years .∈ [ltmp_central[ltmp_central.Year .>= start_year, :Year]])
    comp_years_south = (ref_years .∈ [ltmp_south[ltmp_south.Year .>= start_year, :Year]])

    for ts in axes(res.raw, 1)
        rac = vec(sum(res.raw[ts, :, :], dims=1) .* loc_k_areas') ./ loc_areas
        north_mean_cover[ts] = mean(rac[NORTH_MASK])
        center_mean_cover[ts] = mean(rac[CENTRAL_MASK])
        south_mean_cover[ts] = mean(rac[SOUTH_MASK])
    end

    north_mean_cover = north_mean_cover[comp_years_north]
    center_mean_cover = center_mean_cover[comp_years_center]
    south_mean_cover = south_mean_cover[comp_years_south]

    # Determine performance
    north_perf_series = MAEE_series(north_mean_cover, north_obs)
    central_perf_series = MAEE_series(center_mean_cover, central_obs)
    south_perf_series = MAEE_series(south_mean_cover, south_obs)

    north_perf = temporal_variability(north_perf_series)
    central_perf = temporal_variability(central_perf_series)
    south_perf = temporal_variability(south_perf_series)

    north_t_crp = temporal_correlation(north_perf_series)
    central_t_crp = temporal_correlation(central_perf_series)
    south_t_crp = temporal_correlation(south_perf_series)

    loc_cover = dropdims(sum(res.raw, dims=2), dims=2) .* loc_k_areas' ./ loc_areas'

    class_error_series = class_error(loc_cover)
    reef_error_series = reef_error(loc_cover)

    class_perf = temporal_variability(class_error_series)
    reef_perf = temporal_variability(reef_error_series)

    class_t_crp = temporal_correlation_penalty(series)
    reef_t_crp = temporal_correlation_penalty(series)

    return sum([
        north_perf * north_t_crp,
        central_perf * central_t_crp,
        south_perf * south_t_crp
    ]) / 3 + class_perf * class_t_crp + reef_perf * reef_t_crp
end

base_location_vector = [
    (0.0,  0.95), # habitable cover
    (0.01, 1.0), # taxa weightings
    (0.01, 1.0),
    (0.01, 1.0),
    (0.01, 1.0),
    (0.01, 1.0),
    (0.0,  5.0), # size distribution
    (0.0,  5.0),
    (0.0,  5.0),
    (0.0,  5.0),
    (0.0,  5.0),
]

# Define parameter space to scan over
sample_bounds = repeat(base_location_vector, n_classifications)

best_score_file = "Outputs/best_initial_state_w_coral_$(string(today())).dat"
if !@isdefined(best_init_state) && isfile(best_score_file)
    best_init_state = deserialize(best_score_file)
    @assert all(first.(sample_bounds) .<= best_init_state .<= last.(sample_bounds)) "Initial state is out of bounds"
else
    best_init_state = nothing
end

@info "Using $(Threads.nthreads()-1) threads."
if !@isdefined(best_init_state) || isnothing(best_init_state)
    # Include additional config if using BorgMOEA
    # Method=:borg_moea,
    # FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
    res = bboptimize(
        obj_func;
        SearchRange=sample_bounds,
        MaxSteps=10_000,
        NThreads=Threads.nthreads()-1
    );
elseif !isnothing(best_init_state)
    res = bboptimize(
        obj_func,
        best_init_state;  # provide an initial solution
        SearchRange=sample_bounds,
        MaxSteps=10_000,
        NThreads=Threads.nthreads()-1
    );
end

best_fitness(res)
best_init_state = best_candidate(res)

serialize(best_score_file, best_init_state)
