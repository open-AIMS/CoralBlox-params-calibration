function region_stats(
    region_names::Vector{Symbol},
    region_masks::Vector{BitVector},
    raw_data::Array{Float64,4},
    dom;
    observations::LocationDataStore,
)::NamedTuple
    stats = _region_stats.(Ref(dom), region_masks, Ref(raw_data); observations=observations)
    has_data_mask = .!isnothing.(stats)
    return NamedTuple{Tuple(region_names[has_data_mask])}(Tuple(stats[has_data_mask]))
end

"""
    _region_stats(dom, region_mask, raw_data; observations)

Return aggregated statistics for locations selected by `region_mask`:
- `model_stats`: (3 × n_timesteps) matrix with confint lower, median, upper for model
- `historic_stats`: (3 × n_timesteps) matrix with confint lower, median, upper for observations
- `rmse`, `mae`, `srcc`: pooled metrics across all target locations
"""
function _region_stats(
    dom,
    region_mask::BitVector,
    raw_data::Array{Float64,4};
    observations::LocationDataStore,
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
    @views target_reefs_mask =
        dom.loc_data[region_mask, :].UNIQUE_ID .∈ Ref(observations.ltmp_unique_ids)

    if iszero(target_reefs_mask)
        return nothing
    end

    target_can_reefs_ids = dom.loc_data[region_mask, :][target_reefs_mask, :UNIQUE_ID]
    target_historic = observations.ltmp_coral_cover[
        observations.ltmp_unique_ids .∈ Ref(target_can_reefs_ids), :
    ]

    target_model_data = @view(
        _loc_cover(raw_data, dom)[:, region_mask][:, target_reefs_mask]
    )

    n_locations, n_timesteps = size(target_historic)
    model_stats = Array{Float64}(fill(0.0, 3, n_timesteps))
    historic_stats = Array{Union{Missing,Float64}}(fill(missing, 3, n_timesteps))

    for y in 1:n_timesteps
        model_stats[:, y] .= quantile(@view(target_model_data[y, :]), [0.025 0.5 0.975])'

        target_historic_data = @view(
            target_historic[.!ismissing.(target_historic[:, y]), y]
        )
        isempty(target_historic_data) && continue
        historic_stats[:, y] .= quantile(target_historic_data, [0.025 0.5 0.975])'
    end

    _rmse::Vector{Float64} = zeros(Float64, n_locations)
    _mae::Vector{Float64} = zeros(Float64, n_locations)
    _srcc::Vector{Float64} = zeros(Float64, n_locations)
    for loc_idx in 1:n_locations
        has_data_mask = .!ismissing.(target_historic[loc_idx, :])

        mod_data = target_model_data[:, loc_idx]'[has_data_mask]
        obs_data = Vector{Float64}(target_historic[loc_idx, :][has_data_mask])
        _rmse[loc_idx] = rmse(mod_data, obs_data)
        _mae[loc_idx] = MAE(mod_data, obs_data)
        _srcc[loc_idx] = corspearman(mod_data, obs_data)
    end

    return (
        model_stats=model_stats,
        historic_stats=historic_stats,
        rmse=median(_rmse),
        mae=median(_mae),
        srcc=median(_srcc),
        n_locations=n_locations
    )
end
