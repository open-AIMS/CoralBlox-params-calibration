using ArchGDAL

"""
    calculate_along(lats)

Calculates the normalized position along the GBR from north to south.

# Arguments
- `lats`: Vector of latitudes for each reef

# Returns
Vector of "along" values in [0, 1] where:
- 0 = at the northern end of the GBR
- 1 = at the southern end of the GBR
- 0.5 = midway between north and south

# Formula
along = d_north / (d_north + d_south)

where distances are calculated as latitude differences from the extremes.
"""
function calculate_along(lats)
    # Find the northern and southern ends
    lat_north = maximum(lats)  # Northernmost (highest latitude, less negative)
    lat_south = minimum(lats)  # Southernmost (lowest latitude, more negative)

    # Calculate "along"
    return (lat_north .- lats) ./ (lat_north .- lat_south)
end

## Latitudes
dom_gpkg = VALIDATION_STORE.domain_gpkg

gbr_features = GDF.read("C:/Users/pribeiro/AIMS/Datasets/misc/Great_Barrier_Reef_Features_20_-2075468569012967677.geojson")
mainland_geom = gbr_features[gbr_features.FEAT_NAME.=="Mainland", :geometry][1]

reef_points = ArchGDAL.createpoint.(dom_gpkg.LAT, dom_gpkg.LON)
reef_dists = ArchGDAL.distance.(Ref(mainland_geom), reef_points)

# Step 3
across = reef_dists
along = calculate_along(dom_gpkg.LAT)

# * Corr analysis
valid_idx = [findfirst(dom_gpkg.RME_UNIQUE_ID .== vid) for vid in VALIDATION_STORE.ltmp_unique_ids]
calib_idx = [findfirst(dom_gpkg.RME_UNIQUE_ID .== vid) for vid in CALIBRATION_STORE.ltmp_unique_ids]

across_valid = across[valid_idx]
across_calib = across[calib_idx]

along_valid = along[valid_idx]
along_calib = along[calib_idx]

using Bootstrap
n_boot = 1000

corr_boot = x -> corspearman(x[:, 1], x[:, 2])
bs = BalancedSampling(n_boot)
cil = 0.95

# * Along
Δrmse_c_along_boot = bootstrap(corr_boot, hcat(along_calib, Δrmse_calib), bs)
confint_Δrmse_c_along_boot = confint(Δrmse_c_along_boot, BasicConfInt(cil))

Δrmse_v_along_boot = bootstrap(corr_boot, hcat(along_valid, Δrmse_valid), bs)
confint_Δrmse_v_along_boot = confint(Δrmse_v_along_boot, BasicConfInt(cil))

srcc_c_along_boot = bootstrap(corr_boot, hcat(along_calib, srcc_calib), bs)
confint_srcc_c_along_boot = confint(srcc_c_along_boot, BasicConfInt(cil))

srcc_v_along_boot = bootstrap(corr_boot, hcat(along_valid, srcc_valid), bs)
confint_srcc_v_along_boot = confint(srcc_v_along_boot, BasicConfInt(cil))

# * Across
Δrmse_c_across_boot = bootstrap(corr_boot, hcat(across_calib, Δrmse_calib), bs)
confint_Δrmse_c_across_boot = confint(Δrmse_c_across_boot, BasicConfInt(cil))

Δrmse_v_across_boot = bootstrap(corr_boot, hcat(across_valid, Δrmse_valid), bs)
confint_Δrmse_v_across_boot = confint(Δrmse_v_across_boot, BasicConfInt(cil))

srcc_c_across_boot = bootstrap(corr_boot, hcat(across_calib, srcc_calib), bs)
confint_srcc_c_across_boot = confint(srcc_c_across_boot, BasicConfInt(cil))

srcc_v_across_boot = bootstrap(corr_boot, hcat(across_valid, srcc_valid), bs)
confint_srcc_v_across_boot = confint(srcc_v_across_boot, BasicConfInt(cil))


# * Maximum temperature

max_temps_valid = [maximum(dom.dhw_scens[locs=vid, scenarios=1]) for vid in valid_idx]
max_temps_calib = [maximum(dom.dhw_scens[locs=vid, scenarios=1]) for vid in calib_idx]

Δrmse_c_temp_boot = bootstrap(corr_boot, hcat(max_temps_calib, Δrmse_calib), bs)
confint_Δrmse_c_temp_boot = confint(Δrmse_c_temp_boot, BasicConfInt(cil));

Δrmse_v_temp_boot = bootstrap(corr_boot, hcat(max_temps_valid, Δrmse_valid), bs)
confint_Δrmse_v_temp_boot = confint(Δrmse_v_temp_boot, BasicConfInt(cil));

srcc_c_temp_boot = bootstrap(corr_boot, hcat(max_temps_calib, srcc_calib), bs)
confint_srcc_c_temp_boot = confint(srcc_c_temp_boot, BasicConfInt(cil));

srcc_v_temp_boot = bootstrap(corr_boot, hcat(max_temps_valid, srcc_valid), bs)
confint_srcc_v_temp_boot = confint(srcc_v_temp_boot, BasicConfInt(cil));
