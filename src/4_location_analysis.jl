using CairoMakie, GeoMakie, GraphMakie

include("common/plotting/plotting.jl")
include("1_setup.jl")

corals = ADRIA.to_coral_spec(scens[1, :])

coral_cover = rs_raw.raw .* reshape(site_k_area(dom), (1, 1, 3806))

# Identify temporal regin of interest
temporal_range = 2008:2014

for i in 1:length(limited_locations)
# Location coral params
    limited_loc_pos = i # index of target location in
    location_unique_id = limited_locations[limited_loc_pos]
    domain_loc_pos = findfirst(x->x==location_unique_id, dom.loc_data.UNIQUE_ID)
    ltmp_loc_pos = findfirst(x->!ismissing(x) && x==location_unique_id, ltmp_reef_data.RME_UNIQUE_ID)

    println("---- Plotting Location $(location_unique_id) ----")
    println("Reef Name: $(dom.loc_data.LOC_NAME_L[domain_loc_pos])")
    println("Calibration Location Index: $(limited_loc_pos)")
    println("Ltmp Location Index: $(ltmp_loc_pos)")

    linear_ext = permutedims(reshape(corals.linear_extension, (7, 5)), (2, 1))
    linear_ext[:, 7] .= 0.0
    linear_ext .*= scale_factors[:, limited_loc_pos, 1]

    mb_rate = permutedims(reshape(corals.mb_rate, (7, 5)), (2, 1))
    mb_rate .*= scale_factors[:, limited_loc_pos, 2]

    # Bar plot groupings
    cat = [floor((i - 1) / 5) + 1 for i in 1:35]
    grp = [(i - 1) % 5 + 1 for i in 1:35]

    flt_linear_ext = reshape(linear_ext, (35,))

    save_dir = joinpath(OUT_DIR, "Location_$(location_unique_id)")

    mkpath(save_dir)

    f = plot_coral_param(location_unique_id, "Linear Extension", cat[1:30], grp[1:30], flt_linear_ext[1:30])
    save("$(save_dir)/linexts.png", f)

    flt_mb_rate = reshape(mb_rate, (35,))

    f = plot_coral_param(location_unique_id, "Background Mortality", cat, grp, flt_mb_rate)
    save("$(save_dir)/mbrates.png", f)

    f = plot_taxa_props(location_unique_id, rs_raw.raw[:, :, domain_loc_pos])
    save("$(save_dir)/taxa_props.png", f)

    f = plot_bleaching_mortality(location_unique_id, rs_raw.raw[:, :, domain_loc_pos], rs_raw.bleaching_mortality[:, :, :, domain_loc_pos])
    save("$(save_dir)/bleaching_mortality.png", f)

    f = plot_cyclone_mortality(location_unique_id, rs_raw.raw[:, :, domain_loc_pos], dom.cyclone_mortality_scens[:, domain_loc_pos, :, scens[1, :cyclone_mortality_scenario]].data)
    save("$(save_dir)/cyclone_mortality.png", f)

    f = plot_target_proportions(location_unique_id, rm_ltmp_taxa[:, :, limited_loc_pos])
    save("$(save_dir)/target_composition.png", f)

    xs = dom.loc_data.X_COORD[domain_loc_pos]
    ys = dom.loc_data.Y_COORD[domain_loc_pos]

    #f = ADRIA.viz.map(dom)
    #scatter!([xs], [ys], color=(:black, 0.0), markersize = 20, strokecolor=:black, strokewidth=2, overdraw=true)
    #save("$(save_dir)/loc_map.png", f)

    location_comparison(rs_raw.raw, ltmp_loc_pos, save_dir)
    collect_error_stats(rs_raw.raw, ltmp_loc_pos)

    loc_cov = rs_raw.raw[:, :, domain_loc_pos:domain_loc_pos]

    f = temporal_size_class_proportions(loc_cov)
    save("$(save_dir)/size_class_proportions.png", f)
end
