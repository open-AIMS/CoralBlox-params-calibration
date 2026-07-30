# Generates data/manta_tow_raw_cache.csv: a local cache of ReefMonitoring.get_manta_tow's raw
# output, shared by 03_a_create_cover_manta_csv.jl and generate_clean_growth_intervals.jl so
# neither hits the API independently for the same data. Cached against the full target_loc_ids
# (a superset of what either consumer needs).
#
# Included (guarded on this file's own existence) from base.jl - relies on canonical_gpkg,
# target_loc_ids and DATA_DIR being defined there. Delete data/manta_tow_raw_cache.csv to
# force a fresh pull.

using Dates

@info "No manta_tow_raw_cache.csv found - pulling raw manta tow data from ReefMonitoring"

manta_years_path = joinpath(DATA_DIR, "disturbance_years_manta.csv")
curated_reef_names = isfile(manta_years_path) ?
                      unique(CSV.read(manta_years_path, DataFrame; stringtype=String, comment="#").reef_name) :
                      String[]

spec = ReefMonitoring.get_reef_info()
spec.unique_id = ReefMonitoring.get_unique_id.(spec.latitude, spec.longitude, [canonical_gpkg])
target_reef_names = spec[in.(spec.unique_id, [target_loc_ids]), :aims_reef_name]

reef_names_union = unique(vcat(curated_reef_names, target_reef_names))
println("Pulling manta tow data for ", length(reef_names_union), " reefs...")

manta_tow_cache = DataFrame(reef_name=String[], report_year=Int[], mean=Float64[])
n_errors = 0
for reef_name in reef_names_union
    local manta
    try
        manta = ReefMonitoring.get_manta_tow(reef_name)
    catch e
        global n_errors += 1
        continue
    end
    (isnothing(manta) || isempty(manta)) && continue
    for r in eachrow(manta)
        push!(manta_tow_cache, (reef_name, Int(r.report_year), r.mean))
    end
end
println("Reefs processed: ", length(reef_names_union), " errors: ", n_errors)

out_path = joinpath(DATA_DIR, "manta_tow_raw_cache.csv")
open(out_path, "w") do io
    println(io, "# Raw per-reef, per-year manta tow cover, cached from ReefMonitoring.get_manta_tow.")
    println(io, "# Delete this file and re-run to force a fresh pull.")
    println(io, "# Generated: $(now())")
end
CSV.write(out_path, manta_tow_cache; append=true, writeheader=true)
println("Written to: ", out_path)
