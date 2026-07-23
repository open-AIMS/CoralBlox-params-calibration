using StatsBase

function temporal_correlation(err_series)::Float64
    return cor(1:length(err_series), err_series)
end

function temporal_correlation_penalty(err_series; threshold::Float64=0.3)::Float64
    corr::Float64 = temporal_correlation(err_series)
    return abs(corr) < threshold ? 1.0 : 2 * abs(corr) + 1.0
end

"""
    rmse(modelled, observed)

Calculate Root Mean Squared Error
"""
function rmse(modelled, observed)
    return sqrt(mean((modelled .- observed) .^ 2.0))
end

"""
    bias(modelled, observed)
"""
function bias(modelled, observed)
    return mean(modelled .- observed)
end

"""
Mean Absolute Error
"""
MAE(sim, obs) = mean(abs.(sim .- obs))

"""
Element-wise absolute error, i.e. `MAE` without the final `mean` reduction.
"""
MAE_series(sim, obs) = abs.(sim .- obs)

"""
Mean Absolute Exponential Error.

Assign error that increases exponentially with distance to observed/"true" data.
"""
function MAEE(sim, obs)
    abs_err = abs.(sim .- obs)
    mean(ℯ .^ ((abs_err) .* (1.0 .+ abs_err ./ 1.0)) .- 1.0)
end

function MAEE_series(sim, obs)
    abs_err = abs.(sim .- obs)
    return ℯ .^ ((abs_err) .* (1.0 .+ abs_err / 1.0)) .- 1.0
end

function temporal_variability(x::AbstractVector{<:Real}; w=[0.9, 0.1])
    return mean([mean(x), std(x)], weights(w))
end

function calib_func(sim, obs)
    return temporal_variability(MAEE_series(sim, obs))
end

function _loc_cover(raw_data, dom)
    loc_k_areas = ADRIA.site_k_area(dom)
    loc_areas = ADRIA.loc_area(dom)
    dropdims(sum(raw_data, dims=(2, 3)), dims=(2, 3)) .* loc_k_areas' ./ loc_areas'
end

"""
    collect_error_stats(raw_data, dom; observations)
    collect_error_stats(raw_data, ltmp_loc_idx, dom; observations)

Error stats for all locations or a single location.
"""
function collect_error_stats(
    raw_data::Array{Float64,4},
    dom;
    observations::LocationDataStore,
)::NamedTuple
    n_locations = size(observations.ltmp_coral_cover, 1)

    # Vectors to hold error stats for each location
    rmse_model, rmse_benchmark, mae, pcc, srcc, bias =
        collect.(eachcol(repeat(zeros(Float64, n_locations), 1, n_locations)))

    for loc in 1:n_locations
        error_stats = collect_error_stats(raw_data, loc, dom; observations=observations)

        rmse_model[loc] = error_stats.rmse_model
        rmse_benchmark[loc] = error_stats.rmse_benchmark
        mae[loc] = error_stats.mae
        pcc[loc] = error_stats.pcc
        srcc[loc] = error_stats.srcc
        bias[loc] = error_stats.bias
    end

    error_names = (:rmse_model, :rmse_benchmark, :mae, :pcc, :srcc, :bias)
    return NamedTuple{error_names}((rmse_model, rmse_benchmark, mae, pcc, srcc, bias))
end
function collect_error_stats(
    raw_data::Array{Float64,4},
    ltmp_loc_idx,
    dom;
    observations::LocationDataStore,
)
    loc_cover = _loc_cover(raw_data, dom)

    obs_loc_data = observations.ltmp_coral_cover[ltmp_loc_idx, :]
    not_missing_obs = (!).(ismissing.(obs_loc_data))

    if !any(not_missing_obs)
        return NaN, NaN, NaN, NaN, NaN
    end

    obs_tf = (START_YEAR:END_YEAR)[not_missing_obs]

    domain_idx::Int64 = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    sim_data = loc_cover[:, domain_idx]

    rmse_::Float64 = rmse(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    pcc_::Float64 = cor(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    mae_::Float64 = MAE(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    bias_::Float64 = bias(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    srcc_::Float64 = corspearman(
        sim_data[not_missing_obs], Vector{Float64}(obs_loc_data[not_missing_obs])
    )

    μ_obs = mean(obs_loc_data[not_missing_obs])
    s = length(sim_data[not_missing_obs])
    benchmark_::Float64 = rmse(fill(μ_obs, s), obs_loc_data[not_missing_obs])

    error_names = (:rmse_model, :rmse_benchmark, :mae, :pcc, :srcc, :bias)
    return NamedTuple{error_names}((rmse_, benchmark_, mae_, pcc_, srcc_, bias_))
end

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
    reef_taxa_error(cover; observations)::Float64

Mean compositional mismatch between simulated and observed functional group cover across
LTMP reefs and years. Returns a value in [0, 1] where 0 is a perfect match.

At each (reef, timestep) pair with observations, the Pearson correlation across functional
groups is computed and rescaled from [-1, 1] to [0, 1] via `(1 + r) / 2`, so that:
- r = +1 (perfect match) → 1.0
- r =  0 (no signal)     → 0.5
- r = -1 (inversion)     → 0.0 (worst case, penalised harder than zero correlation)

Per-reef scores are obtained by Fisher z-averaging the per-timestep r values before
rescaling, then the objective is 1 minus the mean score across reefs.
"""
function reef_taxa_error(
    cover;
    observations::LocationDataStore
)
    n_locs = length(observations.composition_to_domain)
    loc_scores = Vector{Float64}(undef, n_locs)

    for (j, idx) in enumerate(observations.composition_to_domain)
        non_missing_mask = vec(all(.!ismissing.(observations.coral_composition[:, :, j]); dims=2))
        observed_ids = eachindex(2008:2022)[non_missing_mask]

        if isempty(observed_ids)
            loc_scores[j] = 0.5
            continue
        end

        r_vals = [
            cor(cover[id, :, idx], Float64.(observations.coral_composition[id, :, j]))
            for id in observed_ids
        ]

        r_mean = tanh(mean(atanh.(r_vals)))
        loc_scores[j] = (1.0 + r_mean) / 2.0
    end

    return 1.0 - mean(loc_scores)
end

"""
    reef_error(cover; observations)::Vector{Float64}

Mean reef-level MAE per year, averaged across all LTMP reefs with observations in that year.

For each LTMP reef location, computes the element-wise MAE between simulated `cover`
and observed coral cover at each observed timestep, then accumulates into a `n_years`
error series. The peak and trough years of each reef's observed time series are weighted
2× to penalise errors at extremes of the dynamic range more heavily. Missing observations
are excluded from both the error sum and the per-year reef count, so each element of the
returned vector is the unweighted mean over whichever reefs had valid data that year.

Returns a `Vector{Float64}` of length `n_years` (one entry per modelled timestep).
"""
function reef_error(
    cover;
    observations::LocationDataStore,
)::Vector{Float64}
    n_years::Int64 = size(cover, 1)
    err_series::Vector{Float64} = zeros(Float64, n_years)
    tmp_err::Vector{Float64} = zeros(Float64, n_years)
    err_counts::Vector{Int64} = zeros(Int64, n_years)
    not_missing::BitVector = BitVector(fill(true, n_years))

    domain_idx::Int64 = -1
    each_obs_ltmp_cover_series = eachrow(observations.ltmp_coral_cover)
    for (row_idx, loc_obs) in enumerate(each_obs_ltmp_cover_series)
        domain_idx = ltmp_cover_idx_to_domain(observations, row_idx)

        not_missing .= (!).(ismissing.(loc_obs))
        trough_idx = argmin(loc_obs[not_missing])
        peak_idx = argmax(loc_obs[not_missing])
        tmp_err[not_missing] .= MAE_series(
            cover[not_missing, domain_idx], loc_obs[not_missing]
        )

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

function rmse_diff(
    rs_raw::Array{Float64,4}, observations::LocationDataStore, dom
)::Vector{Float64}
    n_validation_locs = length(observations.ltmp_cover_to_domain)
    model_rmse = zeros(Float64, n_validation_locs)
    benchmark_rmse = zeros(Float64, n_validation_locs)
    for i in 1:n_validation_locs
        error_stats = collect_error_stats(rs_raw, i, dom; observations=observations)
        model_rmse[i] = error_stats.rmse_model
        benchmark_rmse[i] = error_stats.rmse_benchmark
    end
    return benchmark_rmse .- model_rmse
end

function location_correlation_coefficients(
    raw_data::Array{Float64,4},
    observations::LocationDataStore,
    dom;
    correlation_metric::Symbol=:spearman
)
    n_locs = length(observations.ltmp_cover_to_domain)
    loc_cc = zeros(Float64, n_locs)

    obs_loc_data = observations.ltmp_coral_cover
    not_missing_obs = (!).(ismissing.(obs_loc_data))

    loc_cover = _loc_cover(raw_data, dom)
    domain_idxs = ltmp_cover_idx_to_domain.(Ref(observations), 1:n_locs)
    sim_data = loc_cover[:, domain_idxs]

    correlation_function = if correlation_metric == :spearman
        corspearman
    elseif correlation_metric == :pearson
        cor
    else
        error("Invalid value for `correlation_metric`")
    end

    for i in 1:n_locs
        loc_cc[i] = correlation_function(
            sim_data[not_missing_obs[i, :], i],
            Vector{Float64}(obs_loc_data[i, not_missing_obs[i, :]])
        )
    end

    return loc_cc
end

"""
Use the Fisher transformation to calculate the average correlation coefficient.

Refs:
- https://stats.stackexchange.com/questions/8019/averaging-correlation-values?utm_source=chatgpt.com
- https://www.tandfonline.com/doi/abs/10.1080/00221309809595548
"""
function average_cc(cc_data::Vector{Float64}; w::Vector{Float64}=ones(Float64, length(cc_data)))
    # Apply fisher transformation
    f_data = atanh.(cc_data)

    # Weighted mean
    mean_cc = mean(f_data, weights(w))

    # Apply inverse fisher transformation
    return tanh.(mean_cc)
end
