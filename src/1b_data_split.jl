"""
Randomly split the available data into calibration and validation sets, ensuring there is
both calibration and validation data for every bioregion group.
"""

using Random

# Geopackage containing both the bioregions and the grouped bioregion indices
bioregion_groups_gpkg = GDF.read(BIOREGION_GROUPS_PATH)
# Coral Composition Data
composition_data = open_dataset(COMPOSITION_PATH)


# Load manta tow ltmp reef level data
ltmp_reef_data = GDF.read(LTMP_REEF_DATA_PATH)

# Order year columns in ascending order
ltmp_reef_years = parse.(Int64, names(ltmp_reef_data)[5:end])
ltmp_reef_perm = sortperm(ltmp_reef_years) .+ 4

ltmp_reef_data_names = names(ltmp_reef_data)
ltmp_reef_data_names[5:end] .= ltmp_reef_data_names[ltmp_reef_perm]

# Re-order columns
ltmp_reef_data = select!(ltmp_reef_data, ltmp_reef_data_names...)

# Rescale to be proportions
ltmp_reef_data[:, 5:end] ./= 100

# Remove manta tow observations that did not overlap with a reef polygon or was not
# sufficiently close to a reef polygon
non_missing_mask = (!).(ismissing.(ltmp_reef_data.RME_UNIQUE_ID))
ltmp_reef_data = ltmp_reef_data[non_missing_mask, :]
first_yr_idx = findfirst(x -> x == "2008", names(ltmp_reef_data))

"""
Get all non-unique elements of a list. Found on stack overflow https://stackoverflow.com/questions/54652787/julia-function-to-return-non-unique-elements-of-an-array
"""
function nonunique(x::AbstractArray{T}) where {T}
    xs = sort(x)
    duplicatedvector = T[]
    for i = 2:length(xs)
        if (isequal(xs[i], xs[i-1]) &&
            (length(duplicatedvector) == 0 || !isequal(duplicatedvector[end], xs[i])))

            push!(duplicatedvector, xs[i])
        end
    end

    return duplicatedvector
end

# Use the average coral cover for reef polygons with multiple manta tow observations
non_unique_ltmp_locs = (nonunique(ltmp_reef_data.RME_UNIQUE_ID))
for non_uniq_id in non_unique_ltmp_locs
    # Use the first row found as the row to store the aggregation then delete all others
    first_entry_idx = findfirst(x -> x == non_uniq_id, ltmp_reef_data.RME_UNIQUE_ID)
    non_uniq_mask = ltmp_reef_data.RME_UNIQUE_ID .== non_uniq_id
    non_unique_df = ltmp_reef_data[non_uniq_mask, :]

    # Iterate over each year calculating the mean coral cover
    for col_idx in first_yr_idx:size(ltmp_reef_data, 2)
        if all(ismissing(non_unique_df[:, col_idx]))
            continue
        end
        ltmp_reef_data[first_entry_idx, col_idx] = mean(
            skipmissing(vec(non_unique_df[:, col_idx]))
        )
    end

    non_uniq_mask[first_entry_idx] = false
    to_delete = findall(non_uniq_mask)
    delete!(ltmp_reef_data, to_delete)
end

# Since data has been aggregated between a few locations Reef IDs, GBRMPA IDs and
# geometries are no longer accurate, remove these columns to prevent future errors
select!(ltmp_reef_data, Not(:REEF_ID, :GBRMPA_ID))
first_yr_idx = findfirst(x -> x == "2008", names(ltmp_reef_data))

"""
To minimize the number of index vectors being passed around and copies being created, define
an immutable struct to provide an interface to move between indexing between different
arrays the refer to the same locations.
"""
struct LocationDataStore
    # Data fields
    domain_gpkg::DataFrame
    ltmp_coral_cover::DataFrame
    coral_composition::YAXArray
    # Index Fields
    domain_to_ltmp_cover::Vector{Int64}
    domain_to_composition::Vector{Int64}
    ltmp_cover_to_domain::Vector{Int64}
    composition_to_domain::Vector{Int64}
end

"""
    ltmp_cover_idx_to_domain(loc_data_store::LocationDataStore, ltmp_cover_idx::Int64)

Given an index of a location in the ltmp reef coral cover dataframe, return the index of the
same location in the domain geopackage.
"""
function ltmp_cover_idx_to_domain(loc_data_store::LocationDataStore, ltmp_cover_idx::Int64)
    return loc_data_store.ltmp_cover_to_domain[ltmp_cover_idx]
end

"""
    composition_idx_to_domain(loc_data_store::LocationDataStore, composition_idx::Int64)

Given an index of a location in the coral composition yaxarray, return the index of the
same location in the domain geopackage.
"""
function composition_idx_to_domain(loc_data_store::LocationDataStore, composition_idx::Int64)
    return loc_data_store.composition_to_domain[composition_idx]
end

"""
    domain_idx_to_ltmp_cover(loc_data_store::LocationDataStore, domain_idx::Int64)

Given the index of a location in the domain geopackage, return the index of the same
location in the ltmp reef coral cover dataframe. The majority of domain locations do not
have a ltmp manta tow data.
"""
function domain_idx_to_ltmp_cover(loc_data_store::LocationDataStore, domain_idx::Int64)
    return loc_data_store.domain_to_ltmp_cover[domain_idx]
end

"""
    domain_idx_to_composition(loc_data_store::LocationDataStore, domain_idx::Int64)::Int64

Given the index of a location in the domain geopackage, return the index of the same
location in the coral composition yaxarray. The majority of domain locations do not have
a coral composition data.
"""
function domain_idx_to_composition(loc_data_store::LocationDataStore, domain_idx::Int64)::Int64
    return loc_data_store.domain_to_composition[domain_idx]
end

"""
    sufficient_data(df_row::DataFrameRow)

Ltmp Coral Cover locations with atleast 4 observations that span more then 10 years are
considered to have sufficient data to be used for calibration.
"""
function sufficient_data(df_row::DataFrameRow)
    raw_row = collect(df_row)
    n_obs = count(.!ismissing.(raw_row))
    if n_obs < 4
        return false
    end
    first_obs = findfirst(x -> !ismissing(x), raw_row)
    last_obs = findlast(x -> !ismissing(x), raw_row)
    if last_obs - first_obs < 10
        return false
    end

    return true
end

"""
    split_indices(location_idxs::Vector{Int64}; train_proportion::Float64=0.75)

# Return
- Validation Split
- Calibration Split
"""
function split_indices(location_idxs::Vector{Int64}; calibration_proportion::Float64=0.75)
    n_locations = length(location_idxs)

    n_validation_data::Int64 = Int(floor((1 - calibration_proportion) * n_locations))
    min_validation_data = 2
    if n_validation_data < min_validation_data
        n_validation_data = min_validation_data
    end
    shuffled_idxs = shuffle(location_idxs)

    return shuffled_idxs[1:n_validation_data], shuffled_idxs[n_validation_data+1:end]
end

unique_biogroup_ids = unique(bioregion_groups_gpkg.ASSIGNED_BIOREGION)
sufficient_data_mask = sufficient_data.(eachrow(ltmp_reef_data[:, first_yr_idx:end]))
ltmp_reef_data = ltmp_reef_data[sufficient_data_mask, :]

ltmp_cover_to_domain = [
    findfirst(
        x -> x == uniq_id, dom.loc_data.UNIQUE_ID
    ) for uniq_id in ltmp_reef_data.RME_UNIQUE_ID
]
ltmp_reef_data[!, :BIOGROUP_IDS] .= bioregion_groups_gpkg.ASSIGNED_BIOREGION[ltmp_cover_to_domain]

biogroup_ltmp_idxs = [
    findall(
        ltmp_reef_data.BIOGROUP_IDS .== biogrp_id
    ) for biogrp_id in unique_biogroup_ids
]

biogroup_data_splits = split_indices.(biogroup_ltmp_idxs; calibration_proportion=0.75)

validation_splits = vcat(first.(biogroup_data_splits)...)
calibration_splits = vcat(last.(biogroup_data_splits)...)

validation_df = ltmp_reef_data[validation_splits, :]
calibration_df = ltmp_reef_data[calibration_splits, :]
