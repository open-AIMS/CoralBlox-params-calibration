COLORMAPS = (
    ltmp_disturbances=:seaborn_bright6,
    taxa=:seaborn_colorblind6
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

# Load helpers before everything else
include("./plot_helpers.jl")

include("./location_comparison.jl")
include("./rmse.jl")
include("./progress.jl")
include("./plotting.jl")
