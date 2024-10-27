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

# """
# Asymmetric MAEE which penalizes under-estimates more than over-estimates.
# """
# function MAEE_series(sim, obs)
#     return abs.((1.0 .+ (sim .< obs)) .* ((sim ./ obs) .- 1.0))
# end

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