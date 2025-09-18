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
- `fig_opts` :
    - `:short` : If true, display only RMSE and SCC. Otherwise shows RMSE, μ bnch, PCC, SCC,
    MAEE and BIAS
"""
function _location_err_title(
    raw_data,
    ltmp_loc_idx;
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
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

    err_report_str = if get(fig_opts, :short, true)
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
