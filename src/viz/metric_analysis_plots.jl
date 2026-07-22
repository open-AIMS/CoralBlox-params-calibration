"""
    save_metric_analysis_plots(raw_data, dom, calibration_store, validation_store, out_dir)

Generate and save all metric analysis plots under `out_dir/metrics/`.
"""
function save_metric_analysis_plots(
    raw_data::Array{Float64,4},
    dom,
    calibration_store::LocationDataStore,
    validation_store::LocationDataStore,
    out_dir::String,
)::Nothing
    metrics_save_dir = joinpath(out_dir, "metrics")
    mkpath(metrics_save_dir)

    fig_opts = Dict{Symbol,Any}(:title => "Benchmark RMSE - Model RMSE\nValidation Locations")
    rmse_diff_map = plot_rmse_diff_map(
        raw_data, dom; observations=validation_store, fig_opts=fig_opts
    )
    save(joinpath(metrics_save_dir, "rmse_diff_map.png"), rmse_diff_map)

    fig_opts[:title] = "Pearson Correlation Coefficient\nValidation Locations"
    pearson_coeff_map = plot_correlation_map(
        raw_data, dom; observations=validation_store, metric_type=:pcc, fig_opts=fig_opts
    )
    save(joinpath(metrics_save_dir, "pcc_map.png"), pearson_coeff_map)

    fig_opts[:title] = "Spearman Correlation Coefficient\nValidation Locations"
    scc_map = plot_correlation_map(
        raw_data, dom; observations=validation_store, metric_type=:srcc, fig_opts=fig_opts
    )
    save(joinpath(metrics_save_dir, "scc_map.png"), scc_map)

    rmse_diff_validation = sort(rmse_diff(raw_data, validation_store, dom))
    f_rmse_diff_validation = plot_rmse_scatter(rmse_diff_validation; observation_type="Validation")
    save(joinpath(metrics_save_dir, "rmse_diff_validation.png"), f_rmse_diff_validation)

    rmse_diff_calibration = sort(rmse_diff(raw_data, calibration_store, dom))
    f_rmse_diff_calibration = plot_rmse_scatter(rmse_diff_calibration; observation_type="Calibration")
    save(joinpath(metrics_save_dir, "rmse_diff_calibration.png"), f_rmse_diff_calibration)

    scc_validation = location_correlation_coefficients(raw_data, validation_store, dom)
    f_srcc_validation = plot_metric_scatter(
        scc_validation;
        axis_opts=Dict{Symbol,Any}(
            :title => "Spearman's Correlation Coefficient (SCC)\nValidation reefs",
            :ylabel => "SCC"
        ),
        opts=Dict{Symbol,Any}(:metric_label => "SCC")
    )
    save(joinpath(metrics_save_dir, "scc_validation.png"), f_srcc_validation)

    scc_calibration = location_correlation_coefficients(raw_data, calibration_store, dom)
    f_srcc_calibration = plot_metric_scatter(
        scc_calibration;
        axis_opts=Dict{Symbol,Any}(
            :title => "Spearman's Correlation Coefficient (SCC)\nCalibration reefs",
            :ylabel => "SCC"
        ),
        opts=Dict{Symbol,Any}(:metric_label => "SCC")
    )
    save(joinpath(metrics_save_dir, "scc_calibration.png"), f_srcc_calibration)

    pcc_validation = location_correlation_coefficients(
        raw_data, validation_store, dom; correlation_metric=:pearson
    )
    f_pcc_validation = plot_metric_scatter(
        pcc_validation;
        axis_opts=Dict{Symbol,Any}(
            :title => "Pearson's Correlation Coefficient (PCC)\nValidation reefs",
            :ylabel => "PCC"
        ),
        opts=Dict{Symbol,Any}(:metric_label => "PCC")
    )
    save(joinpath(metrics_save_dir, "pcc_validation.png"), f_pcc_validation)

    pcc_calibration = location_correlation_coefficients(
        raw_data, calibration_store, dom; correlation_metric=:pearson
    )
    f_pcc_calibration = plot_metric_scatter(
        pcc_calibration;
        axis_opts=Dict{Symbol,Any}(
            :title => "Pearson's Correlation Coefficient (PCC)\nCalibration reefs",
            :ylabel => "PCC"
        ),
        opts=Dict{Symbol,Any}(:metric_label => "PCC")
    )
    save(joinpath(metrics_save_dir, "pcc_calibration.png"), f_pcc_calibration)

    fig_m_heatmap = plot_metrics_heatmap(
        raw_data, dom; fig_size=(700, 700), observations=validation_store
    )
    save(joinpath(metrics_save_dir, "metrics_heatmap.png"), fig_m_heatmap)

    return nothing
end
