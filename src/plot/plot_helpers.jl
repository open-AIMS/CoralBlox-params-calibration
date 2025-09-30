"""
    _location_err_title(
        raw_data,
        ltmp_loc_idx;
        fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
        observations::LocationDataStore=COMBINED_STORE,
    )::Makie.RichText

Construct a Rich Test string containing the location unique id, location name and the error
statistics for the location.


# Arguments
- `opts` :
    - `:short_title` : If true, display only RMSE and SCC. Otherwise shows RMSE, μ bnch,
    PCC, SCC, MAEE and BIAS
"""
function _location_err_title(
    raw_data,
    ltmp_loc_idx;
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    observations::LocationDataStore=COMBINED_STORE,
)::Makie.RichText
    domain_idx::Int64 = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    reef_name::String = try
        observations.domain_gpkg.GBR_NAME[domain_idx]
    catch
        observations.domain_gpkg.cluster_id[domain_idx]
    end
    reef_id::String = get_ltmp_loc_unique_id(observations, ltmp_loc_idx)

    rmse_, benchmark_, cc_, maee_, bias_, scc_ = collect_error_stats(
        raw_data, ltmp_loc_idx; observations=observations
    )
    rmse_, benchmark_, cc_, maee_, bias_, scc_ = trunc.(
        [rmse_, benchmark_, cc_, maee_, bias_, scc_], digits=4
    )

    err_report_str = if get(opts, :short_title, true)
        "RMSE: $(rmse_) | SCC: $(scc_)"
    else
        "RMSE: $(rmse_) | μ bnch: $(benchmark_) | " *
        "PCC: $(cc_) | SCC: $(scc_) | MAEE: $(maee_) | BIAS: $(bias_)"
    end

    title_text = rich("$(reef_name)\n$(reef_id)\n", rich(err_report_str, fontsize=15))

    return title_text
end

function fig_coord(idx, n; max_rows=5)
    nrows = min(Int(floor(sqrt(n))), max_rows)
    col = Int(ceil(idx / nrows))
    row = Int(idx - (ceil(idx / nrows) - 1) * nrows)
    return row, col
end

"""
Takes a dictionary `opts`, removes all keys that have one of the substrings in
`except_patterns` and rename all keys that have the pattern `target_pattern` removing this
pattern from the key name.

# Example
```julia
    _opts = Dict{Symbol,Any}(
        :xlabel => "Coral",
        :title => "Coral Cover",
        :xticklabelsvisible_obs_dist => true,
        :yticklabelsvisible_obs_dist => false,
        :xticklabelsvisible_model_dist => false,
        :yticklabelsvisible_model_dist => true
    )

    _new_opts = filter_opts(_opts, obs_dist; except_patterns=["model_dist"])
    # Dict{Symbol,Any}(
    #     :xlabel => "Coral",
    #     :title => "Coral Cover",
    #     :xticklabelsvisible => true,
    #     :yticklabelsvisible => false,
    # )
```
"""
function filter_opts(
    opts::Dict{Symbol,Any}, target_pattern::String; except_patterns::Vector{String}=Vector{String}[]
)::Dict{Symbol,Any}
    regex_pattern = Regex("(" * join(except_patterns, "|") * ")")
    opts_keys = collect(keys(opts))

    # Exclude `except_patterns` keys from `opts`
    new_opt_keys = opts_keys[.!occursin.(regex_pattern, string.(opts_keys))]
    new_opts = filter(((k, v),) -> k ∈ new_opt_keys, opts)

    target_keys = new_opt_keys[occursin.(Ref(target_pattern), string.(new_opt_keys))]
    for key in target_keys
        str_key = string(key)
        new_key = Symbol(strip(replace(str_key, string.(target_pattern) => ""), '_'))
        new_opts[new_key] = pop!(new_opts, key)
    end

    return new_opts
end