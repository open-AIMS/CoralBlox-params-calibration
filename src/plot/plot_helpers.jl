"""
    _location_err_title(raw_data, ltmp_loc_idx; domain_idx=ALL_LTMP_IDXS[ltmp_loc_idx], reef_name=dom.loc_data.GBR_NAME[domain_idx], reef_id=ltmp_reef_data.RME_UNIQUE_ID[ltmp_loc_idx])::Makie.RichText

Construct a Rich Test string containing the location unique id, location name and the error
statistics for the location.
"""
function _location_err_title(
    raw_data,
    ltmp_loc_idx;
    observations::LocationDataStore=COMBINED_STORE,
)::Makie.RichText
    domain_idx::Int64 = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    reef_name::String = observations.domain_gpkg.GBR_NAME[domain_idx]
    reef_id::String = get_ltmp_loc_unique_id(observations, ltmp_loc_idx)

    rmse_, benchmark_, cc_, maee_, bias_ = collect_error_stats(
        raw_data, ltmp_loc_idx; observations=observations
    )
    rmse_, benchmark_, cc_, maee_, bias_ = trunc.(
        [rmse_, benchmark_, cc_, maee_, bias_], digits=4
    )

    err_report_str = "RMSE: $(rmse_) | μ bnch: $(benchmark_) | "
    err_report_str *= "PCC: $(cc_) | MAEE: $(maee_) | BIAS: $(bias_)"

    title_text = rich("$(reef_name)\n$(reef_id)\n", rich(err_report_str, fontsize=15))

    return title_text
end

function fig_coord(idx, n; max_rows=5)
    nrows = min(Int(floor(sqrt(n))), max_rows)
    col = Int(ceil(idx / nrows))
    row = Int(idx - (ceil(idx / nrows) - 1) * nrows)
    return row, col
end
