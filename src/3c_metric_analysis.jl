# Plot Metrics
metrics_save_dir::String = joinpath(OUT_DIR, "metrics")
mkpath(metrics_save_dir)

rmse_diff_map = plot_rmse_diff_map(
    rs_raw.raw;
    observations=VALIDATION_STORE,
    fig_opts=Dict(:title => "Benchmark RMSE - Model RMSE\nValidation Locations")
)
save(joinpath(metrics_save_dir, "rmse_diff_map.png"), rmse_diff_map)

pearson_coeff_map = plot_pcc_map(
    rs_raw.raw;
    observations=VALIDATION_STORE,
    fig_opts=Dict(:title => "Pearson Correlation Coefficient\nValidation Locations")
)
save(joinpath(metrics_save_dir, "pcc_map.png"), pearson_coeff_map)

rmse_diff_validation = sort(rmse_diff(rs_raw.raw, VALIDATION_STORE))
f_rmse_diff_validation = plot_rmse_scatter(rmse_diff_validation; observation_type="Validation")
save(joinpath(metrics_save_dir, "rmse_diff_validation.png"), f_rmse_diff_validation)

rmse_diff_calibration = sort(rmse_diff(rs_raw.raw, CALIBRATION_STORE))
f_rmse_diff_calibration = plot_rmse_scatter(rmse_diff_calibration; observation_type="Calibration")
save(joinpath(metrics_save_dir, "rmse_diff_calibration.png"), f_rmse_diff_calibration)

pcc_validation = pcc_locs(rs_raw.raw, VALIDATION_STORE)
f_pcc_validation = plot_pcc_scatter(pcc_validation; observation_type="Validation")
save(joinpath(metrics_save_dir, "pcc_validation.png"), f_pcc_validation)

pcc_calibration = pcc_locs(rs_raw.raw, CALIBRATION_STORE)
f_pcc_calibration = plot_pcc_scatter(pcc_calibration; observation_type="Calibration")
save(joinpath(metrics_save_dir, "pcc_calibration.png"), f_pcc_calibration)

include("common/common.jl")
fig_m_heatmap = plot_metrics_heatmap(rs_raw.raw; fig_size=(700, 700), observations=VALIDATION_STORE)
save(joinpath(metrics_save_dir, "metrics_heatmap.png"), fig_m_heatmap)
