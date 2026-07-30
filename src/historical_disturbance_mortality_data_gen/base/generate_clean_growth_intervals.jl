# Generates data/clean_growth_intervals.csv: for every CALIBRATION reef (validation reefs
# excluded - they must never feed the growth-rate fit their own mortality is evaluated
# against), every pair of surveyed years (gap 1-6) with NO disturbance of any type reported
# between them, paired with per-functional-group cover at both ends (manta tow total x
# transect composition proportions). A year counts as "surveyed" only if both a manta tow and
# a transect reading exist for it. Feeds 02_estimate_background_growth_rates.jl.
#
# Included (guarded on data/background_growth_rates.csv's existence) from base.jl - relies on
# canonical_gpkg, calibration_loc_ids, CANONICAL_PATH, DATA_DIR and functional_groups being
# defined there. Delete data/background_growth_rates.csv to force regeneration.

using YAXArrays
using Dates

@info "No background_growth_rates.csv found - generating clean_growth_intervals.csv as a prerequisite"

println("Fetching ReefMonitoring reef list...")
spec = ReefMonitoring.get_reef_info()
spec.unique_id = ReefMonitoring.get_unique_id.(spec.latitude, spec.longitude, [canonical_gpkg])

# CALIBRATION-only, not the full target_loc_ids - see README.md, Decisions.
target_spec = spec[in.(spec.unique_id, [calibration_loc_ids]), :]
println("Matched reefs with ReefMonitoring name: ", nrow(target_spec))

# Same frozen coral_composition.nc that 03_b_create_cover_transect_csv.jl uses, not a live
# ReefMonitoring pull - see README.md, Decisions, for why.
println("Loading composition dataset from ", COMPOSITION_PATH)
_raw_composition = open_dataset(COMPOSITION_PATH).mean
composition = YAXArray(
    (
        Dim{:timesteps}(collect(_raw_composition.timesteps.val.data)),
        Dim{:taxa}(replace(
            collect(_raw_composition.taxa.val.data),
            "Corymboses_non_Acropora" => "Corymbose_non_Acropora"
        )),
        Dim{:location}(collect(_raw_composition.location.val.data)),
    ),
    _raw_composition.data,
    _raw_composition.properties
)
composition_locs = Set(collect(composition.location.val.data))

MIN_GAP, MAX_GAP = 1, 6

# Excludes likely unflagged real disturbances (ReefMonitoring's flags can have false
# negatives) - a drop this large in a supposedly clean window is more likely a missed
# disturbance than genuine background decline. Pre-registered, not tuned post-hoc.
SUSPICIOUS_DROP_THRESHOLD = 0.3
MIN_COVER_FOR_DROP_CHECK = 0.005  # avoid div-by-near-zero on reefs with ~0 starting cover

results = DataFrame(
    unique_id=String[], reef_name=String[], year_start=Int[], year_end=Int[], gap_years=Int[],
    n_disturbances_any_type_in_window=Int[], manta_cover_start=Float64[], manta_cover_end=Float64[],
)
for fg in functional_groups
    results[!, "cover_start_$fg"] = Float64[]
    results[!, "cover_end_$fg"] = Float64[]
end

manta_tow_cache = CSV.read(
    joinpath(DATA_DIR, "manta_tow_raw_cache.csv"), DataFrame; stringtype=String, comment="#"
)

n_reef_done = 0
n_errors = 0
n_no_composition = 0
n_no_manta = 0
for row in eachrow(target_spec)
    uid, reef_name = row.unique_id, row.aims_reef_name
    global n_reef_done += 1

    if !(uid in composition_locs)
        global n_no_composition += 1
        continue
    end

    local disturbances
    try
        disturbances = ReefMonitoring.get_disturbances(reef_name)
    catch e
        global n_errors += 1
        continue
    end
    ddates = isnothing(disturbances) || isempty(disturbances) ? Float64[] : parse.(Float64, disturbances.ddate)

    manta = manta_tow_cache[manta_tow_cache.reef_name.==reef_name, [:report_year, :mean]]
    if isempty(manta)
        global n_no_manta += 1
        continue
    end
    manta_by_year = Dict{Int,Float64}()
    for r in eachrow(manta)
        manta_by_year[Int(r.report_year)] = r.mean
    end

    # Bulk-reads the per-reef slice once and indexes in memory - per-scalar disk reads here
    # were found to crash the HDF5 driver after enough reefs.
    loc_cover = composition[location=At(uid)]
    loc_cover_data = Array(loc_cover)
    loc_years = collect(loc_cover.timesteps.val.data)
    loc_taxa = collect(loc_cover.taxa.val.data)
    year_idx = Dict(y => i for (i, y) in enumerate(loc_years))
    taxa_idx = Dict(t => i for (i, t) in enumerate(loc_taxa))

    proportions_by_year = Dict{Int,Vector{Float64}}()
    for yr in 1992:2025
        haskey(year_idx, yr) || continue
        vals = [loc_cover_data[year_idx[yr], taxa_idx[fg]] for fg in functional_groups]
        any(ismissing, vals) && continue
        total = sum(vals)
        total <= 0 && continue
        proportions_by_year[yr] = vals ./ total
    end

    # A year only counts if BOTH manta cover and transect proportions exist for it
    avail_years = sort(collect(intersect(keys(manta_by_year), keys(proportions_by_year))))
    length(avail_years) < 2 && continue

    for i in 1:(length(avail_years) - 1)
        y1 = avail_years[i]
        for j in (i + 1):length(avail_years)
            y2 = avail_years[j]
            gap = y2 - y1
            gap > MAX_GAP && break
            gap < MIN_GAP && continue

            # ddate is a fractional calendar year - compare by calendar year (floor), not the
            # raw fraction, or disturbances after Jan 1st of y2 get silently missed.
            n_dist_in_window = count(d -> y1 <= floor(Int, d) <= y2, ddates)
            n_dist_in_window > 0 && continue  # not clean

            manta_start, manta_end = manta_by_year[y1], manta_by_year[y2]
            if manta_start > MIN_COVER_FOR_DROP_CHECK &&
               (manta_start - manta_end) / manta_start > SUSPICIOUS_DROP_THRESHOLD
                continue  # suspicious drop, likely an unflagged real disturbance
            end

            cover_start = manta_by_year[y1] .* proportions_by_year[y1]
            cover_end = manta_by_year[y2] .* proportions_by_year[y2]

            # NamedTuple keyed by column name, not position - immune to declaration-order bugs
            new_row = merge(
                (unique_id=uid, reef_name=reef_name, year_start=y1, year_end=y2, gap_years=gap,
                    n_disturbances_any_type_in_window=n_dist_in_window,
                    manta_cover_start=manta_by_year[y1], manta_cover_end=manta_by_year[y2]),
                NamedTuple(Symbol("cover_start_$fg") => cover_start[k] for (k, fg) in enumerate(functional_groups)),
                NamedTuple(Symbol("cover_end_$fg") => cover_end[k] for (k, fg) in enumerate(functional_groups)),
            )
            push!(results, new_row)
        end
    end
end

println("Reefs processed: ", n_reef_done, " errors: ", n_errors,
    " no composition data: ", n_no_composition, " no manta data: ", n_no_manta)
println("Clean interval candidates found: ", nrow(results))

sort!(results, [:unique_id, :year_start, :gap_years])

out_path = joinpath(DATA_DIR, "clean_growth_intervals.csv")
open(out_path, "w") do io
    println(io, "# Candidate \"clean\" growth intervals: calibration reefs, consecutive surveyed years with")
    println(io, "# no disturbance of any type between them. cover_start_X/cover_end_X = manta tow total")
    println(io, "# cover x transect composition proportion. Feeds 02_estimate_background_growth_rates.jl.")
    println(io, "# Generated: $(now())")
end
CSV.write(out_path, results; append=true, writeheader=true)
println("Written to: ", out_path)
