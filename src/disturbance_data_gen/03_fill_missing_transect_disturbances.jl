# We use the transect data (benthic composition) to estimate the survival rate for each
# functional group after a disturbance. For each pair disturbance/location without benthic
# data, we use either the average across all locations or across that location's bioregion
# (whichever has the lowest coefficient of variation) to estimate the benthic composition
# before/after the disturbance.

# Load manta and transect data
cover_manta_path = "disturbance_data_gen/data/disturbance_cover_manta.csv"
cover_manta_df = CSV.read(cover_manta_path, DataFrame; stringtype=String, comment="#")

cover_transect_path = "disturbance_data_gen/data/disturbance_cover_transect.csv"
cover_transect_df = CSV.read(cover_transect_path, DataFrame; stringtype=String, comment="#")

# Add `cover_before_<functional_group>` and `cover_after_<functional_group>`
for group in functional_groups
    cover_manta_df[!, "cover_before_$group"] .= 0.0
    cover_manta_df[!, "cover_after_$group"] .= 0.0
end

# For each row in cover_manta_df we need to estimate the benthic composition
# Find manta rows without a corresponding transect row
empty_manta_rows = .!(eachrow(cover_manta_df[!, [:reef_name, :year]]) .∈ [eachrow(cover_transect_df[!, [:reef_name, :year]])])

# Cache Dict to add missing rows to cover_transect_df
cache_cover_row = Dict(pairs(copy(cover_transect_df[1, :])))

for col in Symbol.(vcat(cover_after_col_names, cover_before_col_names, [:year]))
    cache_cover_row[col] = 0.0
end
cache_cover_row[:year_before] = cache_cover_row[:year_after] = cache_cover_row[:reef_id] = 0
cache_cover_row[:reef_name] = cache_cover_row[:cyclone_name] = cache_cover_row[:type] = ""

# For each row in cover_manta_df that doesn't have a corresponding row in cover_transect_df
# add a new corresponding row in cover_transect_df to fill in with cover data
for manta_row in eachrow(cover_manta_df[empty_manta_rows, :])
    cache_cover_row[:cyclone_name] = manta_row.cyclone_name
    cache_cover_row[:type] = manta_row.type
    cache_cover_row[:year] = manta_row.year
    cache_cover_row[:reef_name] = manta_row.reef_name
    cache_cover_row[:reef_id] = manta_row.reef_id
    cache_cover_row[:year_before] = manta_row.year_before
    cache_cover_row[:year_after] = manta_row.year_after
    push!(cover_transect_df, cache_cover_row)
end

# @assert nrow(cover_transect_df) == nrow(cover_manta_df)

sort!(cover_transect_df, [:reef_name, :year])

# Assign bioregions to reefs in cover_transect_df
transect_reef_ids = string.((cover_transect_df.reef_id))
corresponding_indices = [findfirst(==(tid), canonical_gpkg.RME_UNIQUE_ID) for tid in transect_reef_ids]
cover_transect_df.bioregion = Int.(canonical_gpkg[corresponding_indices, :GBRMPA_BIOREGION])

bioregions = sort(unique(cover_transect_df.bioregion))
n_bioregions = length(bioregions)

# Coefficient of variation for each functional group in each bioregion
bioregions_before_cv = zeros(n_bioregions, 5)
bioregions_after_cv = zeros(n_bioregions, 5)

empty_transect_cover = dropdims(sum(Matrix(cover_transect_df[!, cover_before_col_names]), dims=2), dims=2) .== 0.0
# as expected, `empty_manta_rows == empty_transect_cover`

# For each empty row in cover_transect_df we estimate the benthic composition before/after a
# disturbance using either the mean across all reef or the mean inside that reef's bioregion.
# To decide which one to use we compare the Coefficient of Variation (CV) for both sets and
# use the one with smaller CV.

# Coefficient of variation for cover before for each functional group across all reefs
mean_before = mean(Matrix(cover_transect_df[!, cover_before_col_names]), dims=1)
std_before = std(Matrix(cover_transect_df[!, cover_before_col_names]), dims=1)
cv_before = std_before ./ mean_before

mean_after = mean(Matrix(cover_transect_df[!, cover_after_col_names]), dims=1)
std_after = std(Matrix(cover_transect_df[!, cover_after_col_names]), dims=1)
cv_after = std_after ./ mean_after

# For each bioregion fill missing transect cover values with average of bioregion
for (idx, bioregion) in enumerate(bioregions)
    if bioregion == -1
        # Continue unless there's any empty transect cover in this bioregion
        sum(cover_transect_df.bioregion .== bioregion .&& empty_transect_cover) == 0 && continue

        cover_transect_df[cover_transect_df.bioregion.==bioregion.&&empty_transect_cover, cover_before_col_names] .= mean_before
        cover_transect_df[cover_transect_df.bioregion.==bioregion.&&empty_transect_cover, cover_after_col_names] .= mean_after
        continue
    end

    # Rows in cover_transect_df in current bioregion
    target_cover_transect_df = cover_transect_df[cover_transect_df.bioregion.==bioregion, :]

    #
    bio_zero_cover_before = dropdims(sum(Matrix(target_cover_transect_df[:, cover_before_col_names]), dims=2), dims=2) .== 0.0
    bio_zero_cover_after = dropdims(sum(Matrix(target_cover_transect_df[:, cover_after_col_names]), dims=2), dims=2) .== 0.0

    # Locations without cover before should be the same without cover after
    @assert all(bio_zero_cover_before .== bio_zero_cover_after)

    # If there's no data missing in any current bioregion's reef, continue
    if sum(bio_zero_cover_before .&& bio_zero_cover_after) == 0
        continue
    end

    mean_bio_before = mean(Matrix(target_cover_transect_df[.!bio_zero_cover_before, cover_before_col_names]), dims=1)
    std_bio_before = std(Matrix(target_cover_transect_df[.!bio_zero_cover_before, cover_before_col_names]), dims=1)
    bioregions_before_cv[idx, :] = std_bio_before ./ mean_bio_before

    mean_bio_after = mean(Matrix(target_cover_transect_df[.!bio_zero_cover_after, cover_after_col_names]), dims=1)
    std_bio_after = std(Matrix(target_cover_transect_df[.!bio_zero_cover_after, cover_after_col_names]), dims=1)
    bioregions_after_cv[idx, :] = std_bio_after ./ mean_bio_after

    bioregion_cv_before_is_lower = dropdims(bioregions_before_cv[idx, :] .< cv_before', dims=2)
    bioregion_cv_after_is_lower = dropdims(bioregions_after_cv[idx, :] .< cv_after', dims=2)

    replace_cover_before = dropdims(bioregion_cv_before_is_lower' .* mean_bio_before .+ .!bioregion_cv_before_is_lower' .* mean_before, dims=1)
    replace_cover_after = dropdims(bioregion_cv_after_is_lower .* mean_bio_after' .+ .!bioregion_cv_after_is_lower .* mean_after', dims=2)

    cover_transect_df[cover_transect_df.bioregion.==bioregion.&&empty_transect_cover, cover_before_col_names] .= replace_cover_before'
    cover_transect_df[cover_transect_df.bioregion.==bioregion.&&empty_transect_cover, cover_after_col_names] .= replace_cover_after'
end

@assert sum(dropdims(sum(Matrix(cover_transect_df[!, cover_before_col_names]), dims=2), dims=2) .== 0.0) == 0
