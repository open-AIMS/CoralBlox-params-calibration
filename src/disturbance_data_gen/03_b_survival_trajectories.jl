
# Get proportional compositional cover for manta data
reefname_year_transect = cover_transect_df[!, [:reef_name, :year]]

function manta_row_cover(transect_target_cover, manta_cover)
    proportional_transect_cover = transect_target_cover ./ sum(transect_target_cover)
    return (manta_cover .* proportional_transect_cover)'
end

for manta_row in eachrow(cover_manta_df)
    transect_row_mask = eachrow(reefname_year_transect) .== [manta_row[[:reef_name, :year]]]
    transect_target = cover_transect_df[transect_row_mask, :]

    manta_row[cover_before_col_names] .= manta_row_cover(
        Matrix(transect_target[!, cover_before_col_names]), manta_row.cover_before
    )

    sum_before = round(sum(collect(values(manta_row[cover_before_col_names]))), digits=3)
    @assert sum_before == round(manta_row.cover_before, digits=3) || print(sum_before, " ", manta_row.cover_before)

    manta_row[cover_after_col_names] .= manta_row_cover(
        Matrix(transect_target[!, cover_after_col_names]), manta_row.cover_after
    )

    @assert round(sum(collect(values(manta_row[cover_after_col_names]))), digits=3) == round(manta_row.cover_after, digits=3)
end

@assert !any(dropdims(sum(Matrix(cover_manta_df[:, cover_before_col_names]), dims=2) .== 0.0, dims=2))
@assert !any(dropdims(sum(Matrix(cover_manta_df[:, cover_after_col_names]), dims=2) .== 0.0, dims=2))

for group in functional_groups
    survival_rates_df[!, "survival_rate_$group"] .= 0.0
end

survival_rates_cols = "survival_rate_" .* functional_groups

# Calculate survival rates for locations with data
survival_rates_df[:, survival_rates_cols] .= Matrix(cover_manta_df[:, cover_after_col_names]) ./ Matrix(cover_manta_df[:, cover_before_col_names])

# Write
CSV.write("data/survival_rates.csv", survival_rates_df)