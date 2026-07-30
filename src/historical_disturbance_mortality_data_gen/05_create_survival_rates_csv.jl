reefname_year_transect = cover_transect_df[!, [:reef_name, :year]]

function manta_row_cover(transect_target_cover, manta_cover)
    proportional_transect_cover = transect_target_cover ./ sum(transect_target_cover)
    return (manta_cover .* proportional_transect_cover)'
end

transect_proportions(transect_target_cover) = (transect_target_cover ./ sum(transect_target_cover))'

# Caches each record's "before" transect composition proportions (not yet multiplied by any
# cover number) - used below to split the growth-corrected TOTAL cover into per-functional-
# group baselines. See 02_estimate_background_growth_rates.jl for why growth is corrected on
# the total, not per group.
before_proportions = zeros(nrow(cover_manta_df), length(functional_groups))
for (row_idx, manta_row) in enumerate(eachrow(cover_manta_df))
    transect_row_mask = eachrow(reefname_year_transect) .== [manta_row[[:reef_name, :year]]]
    transect_target = cover_transect_df[transect_row_mask, :]

    manta_row[cover_before_col_names] .= manta_row_cover(
        Matrix(transect_target[!, cover_before_col_names]), manta_row.cover_before
    )
    before_proportions[row_idx, :] .= vec(transect_proportions(Matrix(transect_target[!, cover_before_col_names])))

    sum_before = round(sum(collect(values(manta_row[cover_before_col_names]))), digits=3)
    @assert sum_before == round(manta_row.cover_before, digits=3) || print(sum_before, " ", manta_row.cover_before)

    manta_row[cover_after_col_names] .= manta_row_cover(
        Matrix(transect_target[!, cover_after_col_names]), manta_row.cover_after
    )

    sum_after = round(sum(collect(values(manta_row[cover_after_col_names]))), digits=3)
    @assert  sum_after == round(manta_row.cover_after, digits=3)
end

@assert !any(dropdims(sum(Matrix(cover_manta_df[:, cover_before_col_names]), dims=2) .== 0.0, dims=2))
@assert !any(dropdims(sum(Matrix(cover_manta_df[:, cover_after_col_names]), dims=2) .== 0.0, dims=2))

survival_rates_df = copy(cover_manta_df[:, [:reef_name, :reef_id]])

# year_after gave slightly better behavior than year_before when matching disturbance year.
survival_rates_df.dist_year .= cover_manta_df.year_after

for group in functional_groups
    survival_rates_df[!, "survival_rate_$group"] .= 0.0
end

survival_rates_cols = "survival_rate_" .* functional_groups

# Projects TOTAL manta cover forward using the population-level background growth rate,
# compounded over the disturbance's gap (cover_before * (1 + rate)^gap). See
# 02_estimate_background_growth_rates.jl for why growth is modelled on the total, not per
# functional group.
growth_rate_row = only(eachrow(CSV.read(
    joinpath(DATA_DIR, "background_growth_rates.csv"),
    DataFrame; stringtype=String, comment="#"
)))
growth_params = (
    intercept=growth_rate_row.intercept,
    slope_cover_start=growth_rate_row.slope_cover_start,
    slope_gap=growth_rate_row.slope_gap,
)

# Floors the TOTAL compounded factor, not the per-year base - flooring the per-year base at
# 0.0 would compound to exactly 0.0 over any multi-year gap, zeroing the projection and
# turning the survival-rate division into Inf/NaN. Matches the 0.95 mortality ceiling applied
# in 06_build_dist_mortality_datacube.jl.
MIN_TOTAL_GROWTH_FACTOR = 0.05

gaps = cover_manta_df.year_after .- cover_manta_df.year_before
projected_total_cover_before = similar(cover_manta_df.cover_before)
for (row_idx, gap) in enumerate(gaps)
    cover_before_total = cover_manta_df.cover_before[row_idx]
    predicted_rate = growth_params.intercept +
                      growth_params.slope_cover_start * cover_before_total +
                      growth_params.slope_gap * gap
    annual_growth_factor = max(1 + predicted_rate, 0.0)
    total_growth_factor = max(annual_growth_factor^gap, MIN_TOTAL_GROWTH_FACTOR)
    projected_total_cover_before[row_idx] = cover_before_total * total_growth_factor
end

# Split the single growth-corrected TOTAL back into per-functional-group baselines using each
# record's own "before" transect composition weights - the ONLY role transect data plays in
# this projection.
projected_undisturbed_cover = projected_total_cover_before .* before_proportions

# Clamped to [0.05, 1.0], matching the mortality ceiling in
# 06_build_dist_mortality_datacube.jl. Values above 1 are a real, expected residual of a
# population-level correction applied to individual noisy reef trajectories - see README.md,
# Decisions.
raw_survival = Matrix(cover_manta_df[:, cover_after_col_names]) ./ projected_undisturbed_cover

# 0/0 (functional group absent both before and after) can't be fixed by clamping - treated as
# full survival, since there's no coral of that group to have died.
raw_survival[isnan.(raw_survival)] .= 1.0

survival_rates_df[:, survival_rates_cols] .= clamp.(raw_survival, 0.05, 1.0)

CSV.write(joinpath(DATA_DIR, "disturbance_survival_rates.csv"), survival_rates_df)
