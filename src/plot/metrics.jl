using CairoMakie

function plot_rmse_scatter(rmse_diff, observation_type::String)
    _rmse_diff = sort(rmse_diff)

    # Calculate mean and confidence intervals
    mean_rmse_diff = mean(_rmse_diff)
    std_rmse_diff = std(_rmse_diff)
    # ci_lower = mean_rmse_diff - 1.96 * std_rmse_diff / sqrt(length(_rmse_diff))
    # ci_upper = mean_rmse_diff + 1.96 * std_rmse_diff / sqrt(length(_rmse_diff))

    n_greater_than_zero = sum(_rmse_diff .> 0)
    success_rate = round((n_greater_than_zero / length(_rmse_diff)) * 100, digits=2)

    # Create scatter plot
    fig = Figure()
    ax = Axis(fig[1, 1], title="(Benchmark RMSE - Model RMSE) for $(titlecase(observation_type)) Locations\n
    $(success_rate)% of positive diffs", xlabel="Index", ylabel="RMSE Difference")
    scatter!(ax, 1:length(_rmse_diff), _rmse_diff, color=:blue, label="RMSE Diff")

    # Add horizontal line for mean
    hlines!(ax, [mean_rmse_diff], color=:red, linestyle=:dash, label="Mean")

    # Add confidence interval bands
    #band!(ax, collect(1:length(_rmse_diff)),
    #    fill(ci_lower, length(_rmse_diff)), fill(ci_upper, length(_rmse_diff)),
    #    color=(1, 0, 0, 0.8), label="95% CI")

    axislegend(ax)
    fig
end

function plot_pcc_scatter(pcc)
    lower_bound, mean, upper_bound = average_cc(pcc)

    # Create scatter plot
    fig = Figure()
    ax = Axis(fig[1, 1], title="Pearson Correlation Coefficient (PCC) Scatter Plot",
        xlabel="Index", ylabel="PCC")

    scatter!(ax, 1:length(pcc), pcc, color=:blue, label="PCC")

    # Add horizontal line for mean
    hlines!(ax, [mean], color=:red, linestyle=:dash, label="Mean")

    # Add confidence interval bands
    # band!(ax, collect(1:length(pcc)),
    #     fill(lower_bound, length(pcc)), fill(upper_bound, length(pcc)),
    #     color=(0, 0, 1, 0.2), label="Bounds")

    axislegend(ax)
    fig
end
