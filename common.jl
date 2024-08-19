using Revise, Infiltrator
using Serialization

using TOML
using Dates
using CSV, DataFrames, YAXArrays
using StatsBase, Statistics
using BlackBoxOptim

using WGLMakie, GeoMakie, GraphMakie
using ADRIA
using ADRIA: GDF, AG


config = TOML.parsefile("calib_config.toml")
path_configs = config["data_paths"]

# TODO: Move to external config file
reefmod_domain_path = path_configs["reefmod_domain"]
rme_domain_path = path_configs["rme_domain"]
ltmp_shp = path_configs["ltmp_shp"]
canonical_path = path_configs["canonical_path"]
classification_path = path_configs["classification_path"]

location_classification = CSV.read(classification_path, DataFrame)
n_classifications = maximum(location_classification.consecutive_classification)

start_year = 2008
end_year = 2022
n_scenarios = 16

if !@isdefined(OPTIONS)
    uniform_initial_cover = false
    rerun = false

    reefmod_domain = true
    reload_domain = false
    reload_shp = false
    reload_ltmp = false
    reload_canonical = false

    # define OPTIONS to prevent reinitialise on every include
    OPTIONS = true
end

"""
    rmse(modelled, observed)

Calculate Root Mean Squared Error
"""
function rmse(modelled, observed)
    return sqrt(mean((modelled .- observed).^2.0))
end

"""
Mean Absolute Error
"""
MAE(sim, obs) = mean(abs.(sim .- obs))

"""
Mean Absolute Exponential Error.

Assign error that increases exponentially with distance to observed/"true" data.
"""
function MAEE(sim, obs)
    abs_err = abs.(sim .- obs)
    mean(ℯ.^((abs_err ./ obs) .* (1.0 .+ abs_err ./ 1.0)) .- 1.0)
end

function MAEE_series(sim, obs)
    abs_err = abs.(sim .- obs)
    return ℯ.^((abs_err ./ obs) .* (1.0 .+ abs_err / 1.0)) .- 1.0
end

# """
# Asymmetric MAEE which penalizes under-estimates more than over-estimates.
# """
# function MAEE_series(sim, obs)
#     return abs.((1.0 .+ (sim .< obs)) .* ((sim ./ obs) .- 1.0))
# end

function temporal_variability(x::AbstractVector{<:Real}; w=[0.9, 0.1])
    return mean([mean(x), std(x)], weights(w))
end

function calib_func(sim, obs)
    return temporal_variability(MAEE_series(sim, obs))
end


function plot_ltmp(ltmp_n, ltmp_c, ltmp_s)::Nothing
    f = Figure(; size=(900, 1600))
    Axis(
        f[1, 1],
        title="North GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )

    band!(ltmp_n.Year, ltmp_n.lower, ltmp_n.upper, color=(:black, 0.4))
    lines!(ltmp_n.Year, ltmp_n.response, color=:black)
    Axis(
        f[2, 1],
        title="Central GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )


    band!(ltmp_c.Year, ltmp_c.lower, ltmp_c.upper, color=(:black, 0.4))
    lines!(ltmp_c.Year, ltmp_c.response, color=:black)
    Axis(
        f[3, 1],
        title="South GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )

    band!(ltmp_s.Year, ltmp_s.lower, ltmp_s.upper, color=(:black, 0.4))
    lines!(ltmp_s.Year, ltmp_s.response, color=:black)
    save("LTMP_data.png", f)
    return nothing
end

function plot_region(
    f::Figure,
    row::Int64,
    col::Int64,
    title::String,
    obs::DataFrame,
    sim::YAXArray;
    showlegend::Bool = false,
    legend_row::Int64 = 2,
    legend_col::Int64 = 2
)::Axis

    ax = Axis(
        f[row, col],
        title=title,
        xlabel="year",
        ylabel="Relative Absolute Cover",
        width = 400,
        height=400
    )

    mean_agg = dropdims(mean(sim, dims=:sites), dims=:sites)
    confints = ADRIA.analysis.series_confint(mean_agg.data)

    xs::Vector{Float64} = collect(sim.timesteps)
    ADRIA_band = band!(xs, confints[:, 1], confints[:, 3], color=(:red, 0.4))
    ADRIA_line = lines!(xs, confints[:, 2], color=:red)

    xs = obs.Year
    lower::Vector{Float64} = obs.lower
    upper::Vector{Float64} = obs.upper
    response::Vector{Float64} = obs.response

    LTMP_band = band!(xs, lower, upper, color=(:black, 0.4))
    LTMP_line = lines!(xs, response, color=:black)

    if showlegend
        Legend(
            f[legend_row, legend_col],
            [[LTMP_line, LTMP_band], [ADRIA_line, ADRIA_band]],
            ["LTMP", "ADRIAmod"]
        )
    end

    return ax
end

function create_error_statistics(filename::String)::DataFrame
    @info "Computing Error Statistics for time period: $(start_year) - $(end_year)"

    north_ind::Int64 = findfirst(x -> x >= start_year, ltmp_north.Year)
    central_ind::Int64 = findfirst(x -> x >= start_year, ltmp_central.Year)
    south_ind::Int64 = findfirst(x -> x >= start_year, ltmp_south.Year)

    north_end::Int64 = findfirst(x -> x >= end_year, ltmp_north.Year) - 1
    central_end::Int64 = findfirst(x -> x >= end_year, ltmp_central.Year) -1
    south_end::Int64 = findfirst(x -> x >= end_year, ltmp_south.Year) - 1

    north_xs = ltmp_north.Year[north_ind:north_end]
    central_xs = ltmp_central.Year[central_ind:central_end]
    south_xs = ltmp_south.Year[south_ind:south_end]

    north_resp = ltmp_north.response[north_ind:north_end]
    central_resp = ltmp_central.response[central_ind:central_end]
    south_resp = ltmp_south.response[south_ind:south_end]

    s_rac_north = dropdims(mean(north_res[timesteps=At(north_xs)], dims=(:sites, :scenarios)), dims=(:sites, :scenarios))
    s_rac_central = dropdims(mean(central_res[timesteps=At(central_xs)], dims=(:sites, :scenarios)), dims=(:sites, :scenarios))
    s_rac_south = dropdims(mean(south_res[timesteps=At(south_xs)], dims=(:sites, :scenarios)), dims=(:sites, :scenarios))

    rmse_north::Float64 = rmse(s_rac_north, north_resp)
    rmse_central::Float64 = rmse(s_rac_central, central_resp)
    rmse_south::Float64 = rmse(s_rac_south, south_resp)

    @info "RMSE North: $(rmse_north), Central: $(rmse_central), South: $(rmse_south)"

    # Coefficient of Determination
    cc_north::Float64 = cor(s_rac_north, north_resp)
    cc_central::Float64 = cor(s_rac_central, central_resp)
    cc_south::Float64 = cor(s_rac_south, south_resp)

    @info "Correlation Coefficient North: $(cc_north), Central: $(cc_central), South: $(cc_south)"

    err_csv = DataFrame(
        Regions=["North", "Central", "South"],
        RMSE=[rmse_north, rmse_central, rmse_south],
        R=[cc_north, cc_central, cc_south]
    )

    CSV.write(filename, err_csv, writeheader=true)

    return err_csv
end

function plot_residual(
    f::Figure,
    ax_row::Int64,
    ax_col::Int64,
    ltmp_data::DataFrame,
    ADRIA_data::YAXArray,
    title::String
)::Nothing

    Axis(
        f[ax_row, ax_col],
        title=title,
        xlabel="year",
        ylabel="Residuals"
    )

    ltmp_ind::Int64 = findfirst(x -> x >= start_year, ltmp_data.Year)
    ltmp_end::Int64 = findfirst(x -> x >= end_year, ltmp_data.Year) - 1

    ltmp_xs = ltmp_data.Year[ltmp_ind:ltmp_end]
    ltmp_resp = ltmp_data.response[ltmp_ind:ltmp_end]

    adria_mean = dropdims(mean(
        ADRIA_data[timesteps=At(ltmp_xs)], dims=(:sites, :scenarios)
    ), dims=(:sites, :scenarios))

    residuals::Vector{Float64} = ltmp_resp .- adria_mean
    scatter!(ltmp_xs, residuals, color=:black)

    return nothing
end

function plot_residuals(
    filename::String; fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:size => (1600, 1600))
)::Figure

    f = Figure(; fig_opts...)

    plot_residual(
        f,
        1,
        1,
        ltmp_north,
        north_res,
        "North GBR"
    )

    plot_residual(
        f,
        2,
        1,
        ltmp_central,
        central_res,
        "Central GBR"
    )

    plot_residual(
        f,
        3,
        1,
        ltmp_south,
        south_res,
        "South GBR"
    )

    save(filename, f)

    return f
end

function create_compare_plot(
    filename::String; fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:size => (1600, 1600))
)::Figure

    f = Figure(; fig_opts...)

    plot_region(
        f,
        1,
        1,
        "North GBR",
        ltmp_north,
        north_res
    )
    plot_region(
        f,
        1,
        2,
        "Central GBR",
        ltmp_central,
        central_res
    )
    plot_region(
        f,
        2,
        1,
        "South GBR",
        ltmp_south,
        south_res;
        showlegend=true,
        legend_row=2,
        legend_col=2
    )
    resize_to_layout!(f)

    save(filename, f)
    return f
end

# julia> include("path to gist.jl")

# julia> create_compare_plot("plot file name.png")

# julia> create_error_statistics("error statistics filename.csv")

# julia> plot_residuals("residual plot.png")
