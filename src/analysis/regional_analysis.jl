include("../common/perf_metrics.jl")

function region_stats(
    region_names::Vector{Symbol}, region_masks::Vector{BitVector};
    observations::LocationDataStore=COMBINED_STORE
)::NamedTuple
    stats = _region_stats.(Ref(dom), region_masks; observations=observations)
    has_data_mask = .!isnothing.(stats)
    return NamedTuple{Tuple(region_names[has_data_mask])}(Tuple(stats[has_data_mask]))
end

"""
    _region_stats(dom::ADRIA.Domain, region_mask::BitVector;observations::LocationDataStore=COMBINED_STORE)::Union{Nothing,@NamedTuple{model_stats::Array{Union{Missing,Float64}},historic_stats::Array{Union{Missing,Float64}}}}

Return a named tuple with the following aggregated statistics for the locations selected by
the `region_mask`:
    - model_stats: a (3 ⋅ n_timesteps) matrix with confint lower bound, median and confint
upper bound for the modelled predictions
    - historic_stats: a (3 ⋅ n_timesteps) matrix with confint lower bound, median and
confint upper bound for the observed historic data
    - rmse: RMSE for the pooled data across all target locations

"""
function _region_stats(
    dom::ADRIA.Domain,
    region_mask::BitVector;
    observations::LocationDataStore=COMBINED_STORE
)::Union{
    Nothing,
    @NamedTuple{
        model_stats::Array{Float64},
        historic_stats::Array{Union{Missing,Float64}},
        rmse::Float64,
        mae::Float64,
        srcc::Float64,
        n_locations::Int64
    }
}
    target_can_reefs_mask = dom.loc_data[region_mask, :].UNIQUE_ID .∈ Ref(observations.ltmp_unique_ids)

    if iszero(target_can_reefs_mask)
        return nothing
    end

    target_can_reefs_ids = dom.loc_data[region_mask, :][target_can_reefs_mask, :UNIQUE_ID]
    target_historic = observations.ltmp_coral_cover[observations.ltmp_unique_ids.∈Ref(target_can_reefs_ids), :]

    # Select the locations within target region and, from those, the locations within the
    # observation store
    # Also converts rs_raw data from relative to habitable area to relative to total area
    target_model_data = @view(_loc_cover(rs_raw.raw)[:, region_mask][:, target_can_reefs_mask])

    n_locations, n_timesteps = size(target_historic)
    model_stats = Array{Float64}(fill(0.0, 3, n_timesteps))
    historic_stats = Array{Union{Missing,Float64}}(fill(missing, 3, n_timesteps))

    for y in 1:n_timesteps
        # Fill model median and confint
        model_stats[:, y] .= quantile(@view(target_model_data[y, :]), [0.025 0.5 0.975])'

        # Fill historic median and confint
        target_historic_data = @view(target_historic[.!ismissing.(target_historic[:, y]), y])
        isempty(target_historic_data) && continue
        historic_stats[:, y] .= quantile(target_historic_data, [0.025 0.5 0.975])'
    end

    _rmse::Vector{Float64} = zeros(Float64, n_locations)
    _mae::Vector{Float64} = zeros(Float64, n_locations)
    _srcc::Vector{Float64} = zeros(Float64, n_locations)
    _weights::Vector{Float64} = zeros(Float64, n_locations)
    for loc_idx in 1:n_locations
        has_data_mask = .!ismissing.(target_historic[loc_idx, :])
        _weights[loc_idx] = sum(has_data_mask)

        mod_data = target_model_data[:, loc_idx]'[has_data_mask]
        obs_data = Vector{Float64}(target_historic[loc_idx, :][has_data_mask])
        _rmse[loc_idx] = rmse(mod_data, obs_data)
        _mae[loc_idx] = MAE(mod_data, obs_data)
        _srcc[loc_idx] = corspearman(mod_data, obs_data)
    end

    # Mean is weighted by the number of observations (years with data) in each location
    _mean_rmse::Float64 = mean(_rmse, weights(_weights))
    _mean_mae::Float64 = mean(_mae, weights(_weights))
    _mean_srcc::Float64 = average_cc(_srcc; w=_weights)

    return (
        model_stats=model_stats,
        historic_stats=historic_stats,
        rmse=_mean_rmse,
        mae=_mean_mae,
        srcc=_mean_srcc,
        n_locations=n_locations
    )
end
