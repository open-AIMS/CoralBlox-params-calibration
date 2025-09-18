calibration_save_dir::String = joinpath(OUT_DIR, "calibration_locations")
validation_save_dir::String = joinpath(OUT_DIR, "validation_locations")

mkpath(calibration_save_dir)
mkpath(validation_save_dir)

cyc_scens = dom.cyclone_mortality_scens[scenarios=1]
dhw_scens = dom.dhw_scens[scenarios=1]
disturbances_path = "../datasets/ltmp_data/disturbances.nc"
disturbances = open_dataset(disturbances_path).layer

disturbances = open_dataset(disturbances_path).layer

open_dataset("C:/Users/pribeiro/AIMS/Code/CycloneSurvivalRegression/historic_dhw/dhw_scens.nc")

n_calibration_locs = length(CALIBRATION_STORE.ltmp_cover_to_domain)
@showprogress desc = "Plotting calibration locations." for i in 1:n_calibration_locs
    reef_id = CALIBRATION_STORE.ltmp_unique_ids[i]
    f = plot_location_comparison(
        rs_raw.raw, i, dhw_scens, cyc_scens, disturbances;
        fig_opts=Dict{Symbol,Any}(:titlesize => 22),
        hide_ltmp_disturbances=true,
        observations=CALIBRATION_STORE
    )
    save(joinpath(calibration_save_dir, "loc_$(reef_id).png"), f)
end

n_validation_locs = length(VALIDATION_STORE.ltmp_cover_to_domain)
@showprogress desc = "Plotting validation locations." for i in 1:n_validation_locs
    reef_id = VALIDATION_STORE.ltmp_unique_ids[i]
    f = plot_location_comparison(
        rs_raw.raw, i, dhw_scens, cyc_scens, disturbances;
        fig_opts=Dict{Symbol,Any}(:titlesize => 22),
        hide_ltmp_disturbances=true,
        observations=VALIDATION_STORE
    )
    save(joinpath(validation_save_dir, "loc_$(reef_id).png"), f)
end
