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
    return sqrt(mean((modelled .- observed).^2.0))
end

function bias(modelled, observed)
    return mean(modelled .- observed)
end

"""
Mean Absolute Error
"""
MAE(sim, obs) = mean(abs.(sim .- obs))

"""
Mean Absolute Exponential Error.

Assign error that increases exponentially with distance to observed/"true" data.
"""
function MAEE(sim, obs)
    abs_err = abs.(sim .- obs)
    mean(ℯ.^((abs_err) .* (1.0 .+ abs_err ./ 1.0)) .- 1.0)
end

function MAEE_series(sim, obs)
    abs_err = abs.(sim .- obs)
    return ℯ.^((abs_err) .* (1.0 .+ abs_err / 1.0)) .- 1.0
end

function temporal_variability(x::AbstractVector{<:Real}; w=[0.9, 0.1])
    return mean([mean(x), std(x)], weights(w))
end

function calib_func(sim, obs)
    return temporal_variability(MAEE_series(sim, obs))
end

function constant_error_statistics(
    filename,
    stat_func=mean
)::DataFrame
    north_ind::Int64 = findfirst(x -> x >= start_year, ltmp_north.Year)
    central_ind::Int64 = findfirst(x -> x >= start_year, ltmp_central.Year)
    south_ind::Int64 = findfirst(x -> x >= start_year, ltmp_south.Year)

    north_end::Int64 = findfirst(x -> x >= end_year, ltmp_north.Year) - 1
    central_end::Int64 = findfirst(x -> x >= end_year, ltmp_central.Year) -1
    south_end::Int64 = findfirst(x -> x >= end_year, ltmp_south.Year) - 1

    north_xs = ltmp_north.Year[north_ind:north_end]
    central_xs = ltmp_central.Year[central_ind:central_end]
    south_xs = ltmp_south.Year[south_ind:south_end]

    north_resp = ltmp_north.response[north_ind:north_end]
    central_resp = ltmp_central.response[central_ind:central_end]
    south_resp = ltmp_south.response[south_ind:south_end]

    north_stat = repeat([stat_func(north_resp)], length(north_xs))
    central_stat = repeat([stat_func(central_resp)], length(central_xs))
    south_stat = repeat([stat_func(south_resp)], length(south_xs))

    rmse_north::Float64 = rmse(north_stat, north_resp)
    rmse_central::Float64 = rmse(central_stat, central_resp)
    rmse_south::Float64 = rmse(south_stat, south_resp)

    @info "RMSE North: $(rmse_north), Central: $(rmse_central), South: $(rmse_south)"

    # Coefficient of Determination
    cc_north::Float64 = cor(north_stat, north_resp)
    cc_central::Float64 = cor(central_stat, central_resp)
    cc_south::Float64 = cor(south_stat, south_resp)

    @info "Correlation Coefficient North: $(cc_north), Central: $(cc_central), South: $(cc_south)"

    # MAEE
    maee_north::Float64 = MAEE(north_stat, north_resp)
    maee_central::Float64 = MAEE(central_stat, central_resp)
    maee_south::Float64 = MAEE(south_stat, south_resp)

    @info "Mean Absolute Exponential Error North: $(maee_north), Central: $(maee_central), South: $(maee_south)"

    # bias
    bias_north::Float64 = bias(north_stat, north_resp)
    bias_central::Float64 = bias(central_stat, central_resp)
    bias_south::Float64 = bias(south_stat, south_resp)

    @info "Bias North: $(bias_north), Central: $(bias_central), South: $(bias_south)"

    err_csv = DataFrame(
        Regions=["North", "Central", "South"],
        RMSE=[rmse_north, rmse_central, rmse_south],
        R=[cc_north, cc_central, cc_south],
        MAEE=[maee_north, maee_central, maee_south],
        BIAS=[bias_north, bias_central, bias_south]
    )

    CSV.write(filename, err_csv, writeheader=true)

    return err_csv
end

function create_error_statistics(filename::String)::DataFrame
    @info "Computing Error Statistics for time period: $(start_year) - $(end_year)"

    north_ind::Int64 = findfirst(x -> x >= start_year, ltmp_north.Year)
    central_ind::Int64 = findfirst(x -> x >= start_year, ltmp_central.Year)
    south_ind::Int64 = findfirst(x -> x >= start_year, ltmp_south.Year)

    north_end::Int64 = findfirst(x -> x >= end_year, ltmp_north.Year) - 1
    central_end::Int64 = findfirst(x -> x >= end_year, ltmp_central.Year) -1
    south_end::Int64 = findfirst(x -> x >= end_year, ltmp_south.Year) - 1

    north_xs = ltmp_north.Year[north_ind:north_end]
    central_xs = ltmp_central.Year[central_ind:central_end]
    south_xs = ltmp_south.Year[south_ind:south_end]

    north_resp = ltmp_north.response[north_ind:north_end]
    central_resp = ltmp_central.response[central_ind:central_end]
    south_resp = ltmp_south.response[south_ind:south_end]

    s_rac_north = dropdims(mean(north_res[timesteps=At(north_xs)], dims=(:sites, :scenarios)), dims=(:sites, :scenarios))
    s_rac_central = dropdims(mean(central_res[timesteps=At(central_xs)], dims=(:sites, :scenarios)), dims=(:sites, :scenarios))
    s_rac_south = dropdims(mean(south_res[timesteps=At(south_xs)], dims=(:sites, :scenarios)), dims=(:sites, :scenarios))

    rmse_north::Float64 = rmse(s_rac_north, north_resp)
    rmse_central::Float64 = rmse(s_rac_central, central_resp)
    rmse_south::Float64 = rmse(s_rac_south, south_resp)

    @info "RMSE North: $(rmse_north), Central: $(rmse_central), South: $(rmse_south)"

    # Coefficient of Determination
    cc_north::Float64 = cor(s_rac_north, north_resp)
    cc_central::Float64 = cor(s_rac_central, central_resp)
    cc_south::Float64 = cor(s_rac_south, south_resp)

    @info "Correlation Coefficient North: $(cc_north), Central: $(cc_central), South: $(cc_south)"

    # MAEE
    maee_north::Float64 = MAEE(s_rac_north, north_resp)
    maee_central::Float64 = MAEE(s_rac_central, central_resp)
    maee_south::Float64 = MAEE(s_rac_south, south_resp)

    @info "Mean Absolute Exponential Error North: $(maee_north), Central: $(maee_central), South: $(maee_south)"

    # bias
    bias_north::Float64 = bias(s_rac_north, north_resp)
    bias_central::Float64 = bias(s_rac_central, central_resp)
    bias_south::Float64 = bias(s_rac_south, south_resp)

    @info "Bias North: $(bias_north), Central: $(bias_central), South: $(bias_south)"

    err_csv = DataFrame(
        Regions=["North", "Central", "South"],
        RMSE=[rmse_north, rmse_central, rmse_south],
        R=[cc_north, cc_central, cc_south],
        MAEE=[maee_north, maee_central, maee_south],
        BIAS=[bias_north, bias_central, bias_south]
    )

    CSV.write(filename, err_csv, writeheader=true)

    return err_csv
end

function collect_error_stats(
    raw_data,
    ltmp_loc_idx;
    obs_data=all_ltmp_reef,
    obs_idxs=all_ltmp_idxs,
    obs_loc_labels=ltmp_reef_data.RME_UNIQUE_ID,
    loc_k_areas=ADRIA.site_k_area(dom),
    loc_areas=ADRIA.loc_area(dom)
)
    loc_cover = dropdims(sum(raw_data, dims=2), dims=2) .* loc_k_areas' ./ loc_areas'
    if obs_idxs[ltmp_loc_idx] == -1
        return Figure()
    end

    obs_loc_data = obs_data[ltmp_loc_idx, :]
    not_missing_obs = (!).(ismissing.(obs_loc_data))
    obs_tf = (start_year:end_year)[not_missing_obs]

    sim_data = loc_cover[:, obs_idxs[ltmp_loc_idx]]
    reef_id = obs_loc_labels[ltmp_loc_idx]

    rmse_::Float64 = rmse(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    cc_::Float64 = cor(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    maee_::Float64 = MAEE(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    bias_::Float64 = bias(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])

    μ_obs = mean(obs_loc_data[not_missing_obs])
    s = length(sim_data[not_missing_obs])
    benchmark_::Float64 = rmse(fill(μ_obs, s), obs_loc_data[not_missing_obs])

    @info "Location $(reef_id)"
    @info "RMSE:       $(trunc(rmse_, digits=4))"
    @info "Benchmark:  $(trunc(benchmark_, digits=4))"
    @info "Pearsons R: $(trunc(cc_, digits=4))"
    @info "MAEE:       $(trunc(maee_, digits=4))"
    @info "Bias:       $(trunc(bias_, digits=4))"
    return rmse_, benchmark_, cc_, maee_, bias_
end

"""
Calculate the functional group correlation and temporal correlation between aggregated ltmp
data created for reefmod.

Uses the complement of the absolute pearson correlation coefficient such that 0 indicates
a perfect fit, and >= 0 indicates no or negative correlation.

The score indicates the mean correlation.
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
                cor(cover[id, :, idx], rm_ltmp_taxa[id, :, j]) ./ count(non_missing_mask)
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
