using Revise, Infiltrator
using Serialization

using TOML
using Dates
using CSV, DataFrames, YAXArrays
using StatsBase, Statistics
using BlackBoxOptim

using WGLMakie, GeoMakie, GraphMakie
using ADRIA
using ADRIA: GDF, AG, DimensionalData


config = TOML.parsefile("calib_config.toml")
path_configs = config["data_paths"]

# TODO: Move to external config file
reefmod_domain_path = path_configs["reefmod_domain"]
rme_domain_path = path_configs["rme_domain"]
ltmp_shp = path_configs["ltmp_shp"]
canonical_path = path_configs["canonical_path"]
classification_path = path_configs["classification_path"]
manta_tow_class_path = path_configs["manta_tow_path"]
ltmp_reef_data_path = path_configs["ltmp_reef_data"]

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

function taxa_cover_proportions(raw_data)::Figure
    cover = reshape(raw_data, (15, 7, 5, 3806))
    cover = dropdims(sum(cover, dims=4), dims=4)
    cover ./= sum(cover, dims=(2, 3))
    cover = dropdims(sum(cover, dims=2), dims=2)
    cover = permutedims(cover, (2, 1))
    xs = 2008:2022

    f = Figure(; size=(1200, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Proportion",
        title="Functional Group Cover Proportions",
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, cover, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    Legend(f[1, 2], ax, framevisible=false)
    return f
end


function taxa_population_proportions(raw_data)::Figure
    sc_mean_area = reshape(
        permutedims(ADRIA.colony_areas()[2], (2, 1)),
        (1, 35, 1)
    )
    population = raw_data ./ sc_mean_area
    population = reshape(population, (15, 7, 5, 3806))
    population = dropdims(sum(population, dims=4), dims=4)
    population ./= sum(population, dims=(2, 3))
    population = dropdims(sum(population, dims=2), dims=2)
    population = permutedims(population, (2, 1))
    xs = 2008:2022

    f = Figure(; size=(1200, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title="Functional Group Population Proportions",
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, population, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    Legend(f[1, 2], ax, framevisible=false)
    return f
end

function temporal_correlation(series)::Float64
    return cor(1:length(series), series)
end

function temporal_correlation_penalty(series; threshold::Float64=0.3)::Float64
    corr::Float64 = temporal_correlation(series)
    return abs(corr) < threshold ? 1.0 : 2 * abs(corr) + 1.0
end

"""
    temporal_size_class_proportions(raw_data)

Calculate the percentage of coral population occupied by each size class.
"""
function temporal_size_class_proportions(raw_data)::Figure
    sc_mean_area = reshape(
        permutedims(ADRIA.colony_areas()[2], (2, 1)),
        (1, 35, 1)
    )
    population = raw_data ./ sc_mean_area
    population = reshape(population, (15, 7, 5, 3806))
    population = dropdims(sum(population, dims=4), dims=4)
    population ./= sum(population, dims=2)
    population = permutedims(population, (3, 2, 1))

    fg_names = ADRIA.functional_group_names()
    xs = 2008:2022
    col = :PuBuGn_9

    f = Figure(; size=(1600, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[1]),
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, population[1, :, :], color=col, labels="Size Class: " .* string.(1:7))
    Axis(
        f[1, 2];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[2]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[2, :, :], color=col)
    Axis(
        f[2, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[3]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[3, :, :], color=col)
    Axis(
        f[2, 2];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[4]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[4, :, :], color=col)
    Axis(
        f[1, 3];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[5]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[5, :, :], color=col)
    Legend(f[2, 3], ax, framevisible=false)
    colsize!(f.layout, 3, Relative(1/3))
    resize_to_layout!(f)
    return f
end

"""
    rmse(modelled, observed)

Calculate Root Mean Squared Error
"""
function rmse(modelled, observed)
    return sqrt(mean((modelled .- observed).^2.0))
end

function bias(modelled, observed)
    return mean(modelled .- observed)
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
    ADRIA_series = series!(xs, sim[:, :, 1].data[:, :]'; solid_color=(:red, 0.01), linewidth=1, labels=nothing)
    ADRIA_line = lines!(xs, confints[:, 2], color=:red, linewidth=5)

    xs = obs.Year
    lower::Vector{Float64} = obs.lower
    upper::Vector{Float64} = obs.upper
    response::Vector{Float64} = obs.response

    LTMP_band = band!(xs, lower, upper, color=(:black, 0.4))
    LTMP_line = lines!(xs, response, color=:black)

    if showlegend
        Legend(
            f[legend_row, legend_col],
            [[LTMP_line, LTMP_band], [ADRIA_line, ADRIA_series]],
            ["LTMP", "CoralBlox"]
        )
    end

    return ax
end

function constant_error_statistics(
    filename,
    stat_func=mean
)::DataFrame
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

    north_stat = repeat([stat_func(north_resp)], length(north_xs))
    central_stat = repeat([stat_func(central_resp)], length(central_xs))
    south_stat = repeat([stat_func(south_resp)], length(south_xs))

    rmse_north::Float64 = rmse(north_stat, north_resp)
    rmse_central::Float64 = rmse(central_stat, central_resp)
    rmse_south::Float64 = rmse(south_stat, south_resp)

    @info "RMSE North: $(rmse_north), Central: $(rmse_central), South: $(rmse_south)"

    # Coefficient of Determination
    cc_north::Float64 = cor(north_stat, north_resp)
    cc_central::Float64 = cor(central_stat, central_resp)
    cc_south::Float64 = cor(south_stat, south_resp)

    @info "Correlation Coefficient North: $(cc_north), Central: $(cc_central), South: $(cc_south)"

    # MAEE
    maee_north::Float64 = MAEE(north_stat, north_resp)
    maee_central::Float64 = MAEE(central_stat, central_resp)
    maee_south::Float64 = MAEE(south_stat, south_resp)

    @info "Mean Absolute Exponential Error North: $(maee_north), Central: $(maee_central), South: $(maee_south)"

    # bias
    bias_north::Float64 = bias(north_stat, north_resp)
    bias_central::Float64 = bias(central_stat, central_resp)
    bias_south::Float64 = bias(south_stat, south_resp)

    @info "Bias North: $(bias_north), Central: $(bias_central), South: $(bias_south)"

    err_csv = DataFrame(
        Regions=["North", "Central", "South"],
        RMSE=[rmse_north, rmse_central, rmse_south],
        R=[cc_north, cc_central, cc_south],
        MAEE=[maee_north, maee_central, maee_south],
        BIAS=[bias_north, bias_central, bias_south]
    )

    CSV.write(filename, err_csv, writeheader=true)

    return err_csv
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

    # MAEE
    maee_north::Float64 = MAEE(s_rac_north, north_resp)
    maee_central::Float64 = MAEE(s_rac_central, central_resp)
    maee_south::Float64 = MAEE(s_rac_south, south_resp)

    @info "Mean Absolute Exponential Error North: $(maee_north), Central: $(maee_central), South: $(maee_south)"

    # bias
    bias_north::Float64 = bias(s_rac_north, north_resp)
    bias_central::Float64 = bias(s_rac_central, central_resp)
    bias_south::Float64 = bias(s_rac_south, south_resp)

    @info "Bias North: $(bias_north), Central: $(bias_central), South: $(bias_south)"

    err_csv = DataFrame(
        Regions=["North", "Central", "South"],
        RMSE=[rmse_north, rmse_central, rmse_south],
        R=[cc_north, cc_central, cc_south],
        MAEE=[maee_north, maee_central, maee_south],
        BIAS=[bias_north, bias_central, bias_south]
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
