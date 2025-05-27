function region_stats(region_names::Vector{Symbol}, region_masks::Vector{BitVector}; observations::LocationDataStore=COMBINED_STORE)::NamedTuple
    stats = _region_stats.(Ref(dom), region_masks; observations=observations)
    has_data_mask = .!isnothing.(stats)
    return NamedTuple{Tuple(region_names[has_data_mask])}(Tuple(stats[has_data_mask]))
end

"""
    _region_stats(dom::ADRIA.Domain,region_mask::BitVector;observations::LocationDataStore=COMBINED_STORE)::Union{Nothing,@NamedTuple{model_stats::Array{Union{Missing,Float64}},historic_stats::Array{Union{Missing,Float64}}}}

Return a named tuple with model_stats and historic_stats where each is a (3 ⋅ n_timesteps)
Matrix with confint and median values for each timestep
"""
function _region_stats(
    dom::ADRIA.Domain,
    region_mask::BitVector;
    observations::LocationDataStore=COMBINED_STORE
)::Union{
    Nothing,
    @NamedTuple{
        model_stats::Array{Union{Missing,Float64}},
        historic_stats::Array{Union{Missing,Float64}},
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
    #observation store
    target_model = @view(rs_raw.raw[:, :, region_mask][:, :, target_can_reefs_mask])

    # Sum across all functional groups and size classes
    target_model_data = dropdims(sum(target_model, dims=2), dims=2)
    @info target_model_data

    n_locations, n_timesteps = size(target_historic)
    historic_stats = Array{Union{Missing,Float64}}(fill(missing, 3, n_timesteps))
    model_stats = Array{Union{Missing,Float64}}(fill(missing, 3, n_timesteps))

    for y in 1:n_timesteps
        # Fill model median and confint
        @info y, @view(target_model_data[y, :])
        model_stats[:, y] .= quantile(@view(target_model_data[y, :]), [0.025 0.5 0.975])'

        # Fill historic median and confint
        target_historic_data = @view(target_historic[.!ismissing.(target_historic[:, y]), y])
        isempty(target_historic_data) && continue
        historic_stats[:, y] .= quantile(target_historic_data, [0.025 0.5 0.975])'
    end

    return (model_stats=model_stats, historic_stats=historic_stats, n_locations=n_locations)
end
