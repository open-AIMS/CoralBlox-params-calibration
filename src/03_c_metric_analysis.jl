# Plot Metrics
metrics_save_dir::String = joinpath(OUT_DIR, "metrics")
mkpath(metrics_save_dir)

fig_opts = Dict{Symbol,Any}(:title => "")

# Maps
fig_opts[:title] = "Benchmark RMSE - Model RMSE\nValidation Locations"
rmse_diff_map = plot_rmse_diff_map(
    rs_raw.raw; observations=VALIDATION_STORE, fig_opts=fig_opts
)
save(joinpath(metrics_save_dir, "rmse_diff_map.png"), rmse_diff_map)

fig_opts[:title] = "Pearson Correlation Coefficient\nValidation Locations"
pearson_coeff_map = plot_correlation_map(
    rs_raw.raw; observations=VALIDATION_STORE, metric_type=:pcc, fig_opts=fig_opts
)
save(joinpath(metrics_save_dir, "pcc_map.png"), pearson_coeff_map)

fig_opts[:title] = "Spearman Correlation Coefficient\nValidation Locations"
scc_map = plot_correlation_map(
    rs_raw.raw; observations=VALIDATION_STORE, metric_type=:srcc, fig_opts=fig_opts
)
save(joinpath(metrics_save_dir, "scc_map.png"), scc_map)

# Scatter
rmse_diff_validation = sort(rmse_diff(rs_raw.raw, VALIDATION_STORE))
f_rmse_diff_validation = plot_rmse_scatter(rmse_diff_validation; observation_type="Validation")
save(joinpath(metrics_save_dir, "rmse_diff_validation.png"), f_rmse_diff_validation)

rmse_diff_calibration = sort(rmse_diff(rs_raw.raw, CALIBRATION_STORE))
f_rmse_diff_calibration = plot_rmse_scatter(rmse_diff_calibration; observation_type="Calibration")
save(joinpath(metrics_save_dir, "rmse_diff_calibration.png"), f_rmse_diff_calibration)

scc_validation = location_correlation_coefficients(rs_raw.raw, VALIDATION_STORE)
f_srcc_validation = plot_metric_scatter(
    scc_validation;
    axis_opts=Dict{Symbol,Any}(
        :title => "Spearman's Correlation Coefficient (SCC)\nValidation reefs",
        :ylabel => "SCC"
    ),
    opts=Dict{Symbol,Any}(:metric_label => "SCC")
)
save(joinpath(metrics_save_dir, "scc_validation.png"), f_srcc_validation)

scc_calibration = location_correlation_coefficients(rs_raw.raw, CALIBRATION_STORE)
f_srcc_calibration = plot_metric_scatter(
    scc_calibration;
    axis_opts=Dict{Symbol,Any}(
        :title => "Spearman's Correlation Coefficient (SCC)\nCalibration reefs",
        :ylabel => "SCC"
    ),
    opts=Dict{Symbol,Any}(:metric_label => "SCC")
)
save(joinpath(metrics_save_dir, "scc_calibration.png"), f_srcc_calibration)

pcc_validation = location_correlation_coefficients(rs_raw.raw, VALIDATION_STORE; correlation_metric=:pearson)
f_pcc_validation = plot_metric_scatter(
    pcc_validation;
    axis_opts=Dict{Symbol,Any}(
        :title => "Pearson's Correlation Coefficient (PCC)\nValidation reefs",
        :ylabel => "PCC"
    ),
    opts=Dict{Symbol,Any}(:metric_label => "PCC")
)
save(joinpath(metrics_save_dir, "pcc_validation.png"), f_pcc_validation)

pcc_calibration = location_correlation_coefficients(rs_raw.raw, CALIBRATION_STORE; correlation_metric=:pearson)
f_pcc_calibration = plot_metric_scatter(
    pcc_calibration;
    axis_opts=Dict{Symbol,Any}(
        :title => "Pearson's Correlation Coefficient (PCC)\nCalibration reefs",
        :ylabel => "PCC"
    ),
    opts=Dict{Symbol,Any}(:metric_label => "PCC")
)
save(joinpath(metrics_save_dir, "pcc_calibration.png"), f_pcc_calibration)

# Heatmaps
fig_m_heatmap = plot_metrics_heatmap(rs_raw.raw; fig_size=(700, 700), observations=VALIDATION_STORE)
save(joinpath(metrics_save_dir, "metrics_heatmap.png"), fig_m_heatmap)
