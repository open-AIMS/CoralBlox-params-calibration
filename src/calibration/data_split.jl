struct CalibrationData
    calibration_store::LocationDataStore
    validation_store::LocationDataStore
    combined_store::LocationDataStore
end

function _nonunique(x::AbstractArray{T}) where {T}
    xs = sort(x)
    dups = T[]
    for i in 2:length(xs)
        if isequal(xs[i], xs[i-1]) && (isempty(dups) || !isequal(dups[end], xs[i]))
            push!(dups, xs[i])
        end
    end
    return dups
end

function _sufficient_data(df_row::DataFrameRow)::Bool
    raw_row = collect(df_row)
    n_obs = count(.!ismissing.(raw_row))
    n_obs < 4 && return false
    first_obs = findfirst(x -> !ismissing(x), raw_row)
    last_obs = findlast(x -> !ismissing(x), raw_row)
    return last_obs - first_obs >= 10
end

function _split_indices(location_idxs::Vector{Int64}; calibration_proportion::Float64=0.75)
    n_validation = max(2, Int(floor((1 - calibration_proportion) * length(location_idxs))))
    shuffled = shuffle(location_idxs)
    return shuffled[1:n_validation], shuffled[n_validation+1:end]
end

"""
    build_calibration_data(dom, ltmp_reef_data_path, composition_path; out_dir, rng_seed)

Load and split LTMP manta-tow and coral composition data into calibration, validation, and
combined `LocationDataStore`s. Writes `calibration_split.csv` to `out_dir` when provided.
"""
function build_calibration_data(
    dom,
    ltmp_reef_data_path::String,
    composition_path::String;
    out_dir::Union{String,Nothing}=nothing,
    rng_seed::Int=1
)::CalibrationData
    composition_data = open_dataset(composition_path)

    ltmp_reef_data = GDF.read(ltmp_reef_data_path)

    # Sort year columns ascending
    ltmp_reef_data_names = names(ltmp_reef_data)
    start_year_column = findfirst(ltmp_reef_data_names .== "1993")
    year_columns = start_year_column:lastindex(ltmp_reef_data_names)
    ltmp_reef_years = parse.(Int64, ltmp_reef_data_names[year_columns])
    ltmp_reef_perm = sortperm(ltmp_reef_years) .+ (start_year_column - 1)
    ltmp_reef_data_names[year_columns] .= ltmp_reef_data_names[ltmp_reef_perm]
    ltmp_reef_data = select!(ltmp_reef_data, ltmp_reef_data_names...)
    ltmp_reef_data[:, year_columns] ./= 100

    # Remove missing / duplicate locations
    ltmp_reef_data = ltmp_reef_data[(!).(ismissing.(ltmp_reef_data.RME_UNIQUE_ID)), :]
    first_yr_idx = findfirst(x -> x == "2008", names(ltmp_reef_data))

    for non_uniq_id in _nonunique(ltmp_reef_data.RME_UNIQUE_ID)
        first_entry_idx = findfirst(x -> x == non_uniq_id, ltmp_reef_data.RME_UNIQUE_ID)
        non_uniq_mask = ltmp_reef_data.RME_UNIQUE_ID .== non_uniq_id
        non_unique_df = ltmp_reef_data[non_uniq_mask, :]
        for col_idx in first_yr_idx:size(ltmp_reef_data, 2)
            all(ismissing.(non_unique_df[:, col_idx])) && continue
            ltmp_reef_data[first_entry_idx, col_idx] = mean(
                skipmissing(vec(non_unique_df[:, col_idx]))
            )
        end
        non_uniq_mask[first_entry_idx] = false
        delete!(ltmp_reef_data, findall(non_uniq_mask))
    end

    select!(ltmp_reef_data, Not(:REEF_ID, :GBRMPA_ID, :geometry))
    select!(ltmp_reef_data, [:RME_UNIQUE_ID, Symbol.(START_YEAR:END_YEAR)...])
    first_yr_idx = findfirst(x -> x == "2008", names(ltmp_reef_data))

    sufficient_data_mask = _sufficient_data.(eachrow(ltmp_reef_data[:, first_yr_idx:end]))
    ltmp_reef_data = ltmp_reef_data[sufficient_data_mask, :]

    # Map LTMP locations to domain indices and assign biogroup IDs
    unique_biogroup_ids = unique(dom.loc_data.CB_CALIB_GROUPS)
    ltmp_cover_to_domain = [
        findfirst(x -> x == id, dom.loc_data.UNIQUE_ID)
        for id in ltmp_reef_data.RME_UNIQUE_ID
    ]
    ltmp_reef_data[!, :BIOGROUP_IDS] .= dom.loc_data.CB_CALIB_GROUPS[ltmp_cover_to_domain]

    biogroup_ltmp_idxs = [
        findall(ltmp_reef_data.BIOGROUP_IDS .== bg) for bg in unique_biogroup_ids
    ]

    Random.seed!(rng_seed)
    biogroup_data_splits = _split_indices.(biogroup_ltmp_idxs; calibration_proportion=0.75)
    validation_splits = vcat(first.(biogroup_data_splits)...)
    calibration_splits = vcat(last.(biogroup_data_splits)...)

    select!(ltmp_reef_data, Not(:BIOGROUP_IDS))
    validation_df = ltmp_reef_data[validation_splits, :]
    calibration_df = ltmp_reef_data[calibration_splits, :]

    cal_ltmp_to_domain = [
        findfirst(dom.loc_data.UNIQUE_ID .== id) for id in calibration_df.RME_UNIQUE_ID
    ]
    val_ltmp_to_domain = [
        findfirst(dom.loc_data.UNIQUE_ID .== id) for id in validation_df.RME_UNIQUE_ID
    ]
    all_ltmp_to_domain = [
        findfirst(dom.loc_data.UNIQUE_ID .== id) for id in ltmp_reef_data.RME_UNIQUE_ID
    ]

    calibration_comp_ids = [
        id for id in composition_data.location if !(id in validation_df.RME_UNIQUE_ID)
    ]
    validation_comp_ids = [
        id for id in composition_data.location if id in validation_df.RME_UNIQUE_ID
    ]
    calibration_mask = (!).(in.(composition_data.location, Ref(validation_df.RME_UNIQUE_ID)))
    validation_mask = in.(composition_data.location, Ref(validation_df.RME_UNIQUE_ID))

    cal_composition_to_domain = [
        findfirst(dom.loc_data.UNIQUE_ID .== id) for id in calibration_comp_ids
    ]
    val_composition_to_domain = [
        findfirst(dom.loc_data.UNIQUE_ID .== id) for id in validation_comp_ids
    ]
    all_composition_to_domain = [
        findfirst(dom.loc_data.UNIQUE_ID .== id) for id in composition_data.location
    ]

    if !isnothing(out_dir)
        calibration_unique_ids = dom.loc_data.UNIQUE_ID[cal_ltmp_to_domain]
        validation_unique_ids = dom.loc_data.UNIQUE_ID[val_ltmp_to_domain]
        calib_valid_split = DataFrame(
            :UNIQUE_IDS => vcat(calibration_unique_ids, validation_unique_ids),
            :BIOGROUP => vcat(
                dom.loc_data.CB_CALIB_GROUPS[cal_ltmp_to_domain],
                dom.loc_data.CB_CALIB_GROUPS[val_ltmp_to_domain]
            ),
            :USAGE => vcat(
                fill("calibration", length(cal_ltmp_to_domain)),
                fill("validation", length(val_ltmp_to_domain))
            )
        )
        split_path = joinpath(out_dir, "calibration_split.csv")
        CSV.write(split_path, calib_valid_split)
        @info "Saving calibration/validation split to $split_path"
    end

    calibration_store = LocationDataStore(
        dom.loc_data,
        calibration_df.RME_UNIQUE_ID,
        Matrix(calibration_df[:, 2:end]),
        composition_data.mean[location=calibration_mask, timestep=At(START_YEAR:END_YEAR)].data[:, :, :],
        cal_ltmp_to_domain,
        cal_composition_to_domain
    )
    validation_store = LocationDataStore(
        dom.loc_data,
        validation_df.RME_UNIQUE_ID,
        Matrix(validation_df[:, 2:end]),
        composition_data.mean[location=validation_mask, timestep=At(START_YEAR:END_YEAR)].data[:, :, :],
        val_ltmp_to_domain,
        val_composition_to_domain
    )
    combined_store = LocationDataStore(
        dom.loc_data,
        ltmp_reef_data.RME_UNIQUE_ID,
        Matrix(ltmp_reef_data[:, 2:end]),
        composition_data.mean[timestep=At(START_YEAR:END_YEAR)].data[:, :, :],
        all_ltmp_to_domain,
        all_composition_to_domain
    )

    return CalibrationData(calibration_store, validation_store, combined_store)
end
