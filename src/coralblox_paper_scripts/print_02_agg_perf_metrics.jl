using Bootstrap

# * Bootstrap setup
n_boot = 1000
bs = BalancedSampling(n_boot)
cil = 0.95

# * Bootstrapped performance metrics aggregates

# ** Calib

stats_calib = collect_error_stats(rs_raw.raw; observations=CALIBRATION_STORE)

rmse_model_c = stats_calib.rmse_model
rmse_benchmark_c = stats_calib.rmse_benchmark
Δrmse_c = rmse_benchmark_c .- rmse_model_c
srcc_model_c = stats_calib.srcc

mean_rmse_c_boot = bootstrap(mean, rmse_model_c, bs)
confint_rmse_c_boot = confint(mean_rmse_c_boot, BasicConfInt(cil))

mean_Δrmse_c_boot = bootstrap(mean, Δrmse_c, bs)
confint_Δrmse_c_boot = confint(mean_Δrmse_c_boot, BasicConfInt(cil))

mean_srcc_c_boot = bootstrap(mean, srcc_model_c, bs)
confint_srcc_c_boot = confint(mean_srcc_c_boot, BasicConfInt(cil))

# ** Valid
stats_valid = collect_error_stats(rs_raw.raw; observations=VALIDATION_STORE)
rmse_model_v = stats_valid.rmse_model
rmse_benchmark_v = stats_valid.rmse_benchmark
Δrmse_v = rmse_benchmark_v .- rmse_model_v
srcc_model_v = stats_valid.srcc

mean_rmse_v_boot = bootstrap(mean, rmse_model_v, bs)
confint_rmse_v_boot = confint(mean_rmse_v_boot, BasicConfInt(cil))

mean_Δrmse_v_boot = bootstrap(mean, Δrmse_v, bs)
confint_Δrmse_v_boot = confint(mean_Δrmse_v_boot, BasicConfInt(cil))

mean_srcc_v_boot = bootstrap(mean, srcc_model_v, bs)
confint_srcc_v_boot = confint(mean_srcc_v_boot, BasicConfInt(cil))