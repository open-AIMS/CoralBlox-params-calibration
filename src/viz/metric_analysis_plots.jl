"""
    save_metric_analysis_plots(raw_data, dom, calibration_store, test_store, out_dir)

Generate and save all metric analysis plots under `out_dir/metrics/`.
"""
function save_metric_analysis_plots(
    raw_data::Array{Float64,4},
    dom,
    calibration_store::LocationDataStore,
    test_store::LocationDataStore,
    out_dir::String,
)::Nothing
    metrics_save_dir = joinpath(out_dir, "metrics")
    mkpath(metrics_save_dir)

    fig_opts = Dict{Symbol,Any}(:title => "Benchmark RMSE - Model RMSE\nTest Locations")
    rmse_diff_map = plot_rmse_diff_map(
        raw_data, dom; observations=test_store, fig_opts=fig_opts
    )
    save(joinpath(metrics_save_dir, "rmse_diff_map.png"), rmse_diff_map)

    fig_opts[:title] = "Pearson Correlation Coefficient\nTest Locations"
    pearson_coeff_map = plot_correlation_map(
        raw_data, dom; observations=test_store, metric_type=:pcc, fig_opts=fig_opts
    )
    save(joinpath(metrics_save_dir, "pcc_map.png"), pearson_coeff_map)

    fig_opts[:title] = "Spearman Correlation Coefficient\nTest Locations"
    scc_map = plot_correlation_map(
        raw_data, dom; observations=test_store, metric_type=:srcc, fig_opts=fig_opts
    )
    save(joinpath(metrics_save_dir, "scc_map.png"), scc_map)

    rmse_diff_stats_test = rmse_diff_stats(raw_data, test_store, dom)
    f_rmse_diff_test = plot_rmse_scatter(
        rmse_diff_stats_test; observation_type="Test"
    )
    save(joinpath(metrics_save_dir, "rmse_diff_test.png"), f_rmse_diff_test)

    rmse_diff_stats_calibration = rmse_diff_stats(raw_data, calibration_store, dom)
    f_rmse_diff_calibration = plot_rmse_scatter(
        rmse_diff_stats_calibration; observation_type="Calibration"
    )
    save(joinpath(metrics_save_dir, "rmse_diff_calibration.png"), f_rmse_diff_calibration)

    nse_stats_test = nse_stats(raw_data, test_store, dom)
    f_nse_test = plot_nse_scatter(nse_stats_test; observation_type="Test")
    save(joinpath(metrics_save_dir, "nse_test.png"), f_nse_test)

    nse_stats_calibration = nse_stats(raw_data, calibration_store, dom)
    f_nse_calibration = plot_nse_scatter(nse_stats_calibration; observation_type="Calibration")
    save(joinpath(metrics_save_dir, "nse_calibration.png"), f_nse_calibration)

    srcc_stats_test = correlation_stats(
        raw_data, test_store, dom; correlation_metric=:spearman
    )
    f_srcc_test = plot_correlation_scatter(
        srcc_stats_test; metric_label="SRCC", observation_type="Test"
    )
    save(joinpath(metrics_save_dir, "scc_test.png"), f_srcc_test)

    srcc_stats_calibration = correlation_stats(
        raw_data, calibration_store, dom; correlation_metric=:spearman
    )
    f_srcc_calibration = plot_correlation_scatter(
        srcc_stats_calibration; metric_label="SRCC", observation_type="Calibration"
    )
    save(joinpath(metrics_save_dir, "scc_calibration.png"), f_srcc_calibration)

    pcc_stats_test = correlation_stats(
        raw_data, test_store, dom; correlation_metric=:pearson
    )
    f_pcc_test = plot_correlation_scatter(
        pcc_stats_test; metric_label="PCC", observation_type="Test"
    )
    save(joinpath(metrics_save_dir, "pcc_test.png"), f_pcc_test)

    pcc_stats_calibration = correlation_stats(
        raw_data, calibration_store, dom; correlation_metric=:pearson
    )
    f_pcc_calibration = plot_correlation_scatter(
        pcc_stats_calibration; metric_label="PCC", observation_type="Calibration"
    )
    save(joinpath(metrics_save_dir, "pcc_calibration.png"), f_pcc_calibration)

    fig_m_heatmap = plot_metrics_heatmap(
        raw_data, dom; fig_size=(700, 700), observations=test_store
    )
    save(joinpath(metrics_save_dir, "metrics_heatmap.png"), fig_m_heatmap)

    return nothing
end
