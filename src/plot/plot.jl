using CairoMakie
using ColorSchemes

COLORMAPS = (
    ltmp_disturbances=:seaborn_bright6,
    taxa=:seaborn_colorblind6
)

COLORS = (
    model_vs_obs_color_model="#1E88E5",
    model_vs_obs_color_obs="#D81B60",
    model_dist_color_dhw="#eb4034",
    model_dist_color_dist="#31004D",
    model_line_color="#1E88E5",
    model_band_color=("#1E88E5", 0.4),
    historic_line_color="#D81B60",
    historic_band_color=("#D81B60", 0.3),
    NamedTuple(zip(
        Symbol.(ADRIA.functional_group_names()),
        ColorSchemes.seaborn_colorblind6[1:5],
    )
    )...
)

FONT_SIZES = (
    xlabel=20,
    ylabel=20,
    title=20
)

BASE_AXIS_OPTS::Dict = Dict(
    :xlabelsize => FONT_SIZES.xlabel,
    :ylabelsize => FONT_SIZES.ylabel,
    :xticks => collect(START_YEAR:END_YEAR),
    :xticklabelrotation => pi / 6,
    :limits => ((START_YEAR - 0.5, END_YEAR + 0.5), nothing),
)

FIG_SIZE::Dict = Dict(
    :map => (500, 600),
)

GBRMPA_MAINLAND_GPKG = GDF.read(GBRMPA_MAINLAND_PATH)

# * Units relative to 1 CSS px
inch = 96
pt = 4 / 3
cm = inch / 2.54

# Load helpers before everything else
include("./plot_helpers.jl")

include("./location_comparison.jl")
include("./regional_comparison.jl")
include("./metrics.jl")
include("./progress.jl")
include("./plotting.jl")

@info("Plots loaded")
