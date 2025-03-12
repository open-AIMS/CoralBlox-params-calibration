# ADRIA-CoralBlox calibration

## Setup

Instantiate environment to install required packages.

```julia
]instantiate
```

It is assumed that all scripts are run inside the `src` directory.

Add the required version of ADRIA using:
`add https://www.github.com/open-AIMS/ADRIA.jl#takuya/calib`

The latest version of CoralBlox is also required:
`add https://www.github.com/open-AIMS/CoralBlox.jl#main`

Create a `calib_config.toml` file with the following entries:

```toml
[Domains]
reefmod_domain = "<path to ReefMod dataset>"
rme_domain = "<path to RME dataset>"
historic_cyclone_mortality = "<path to historic cyclone mortality dataset>"

[Geospatial]
canonical_path = "<path to canonical gpkg>"
ltmp_shp = "../spatial_data/gbr_3Zone 2.shp"
classification_path = "../spatial_data/location_classification_MPA.csv"

[Observations]
manta_tow_path = "..\\ltmp_data\\manta_tow_mean_std.nc"
ltmp_reef_data = "..\\ltmp_data\\manta_tow_data_reef_lvl.gpkg"
composition_netcdf = "..\\ltmp_data\\coral_composition.nc"
ltmp_modelled_obs = "..\\ltmp_data\\modelled_brms.beta.ry.disp.csv"

[Initialisation]
init_cover_filepath = "spatial_data\\init_cover.dat"
init_guess_filepath = "<path to initial guess.dat>" # optional

[Outputs]
out_dir = "..\\Outputs\\test_dir"
result_filename = "results.dat" # relative to out_dir
```

## Datasets

### RME Data

#### [rme_cover_20_reps_2008_2022.nc](datasets/rme_data/rme_cover_20_reps_2008_2022.nc)

Results for 20 ReefModEngine repetitions. All repetitions use the same historic DHW scenario but were initialized with distinct initial coral cover values. These were taken from ReefMod-GBR CMIP6 counterfactuals (Apr. 2024), GCM CNRM-ESM2-1 and DHW SSP1-1.9. The original data was sliced to the timeframe from 2008 to 2022. Link: [https://data.mds.gbrrestoration.org/dataset/102.100.100/653201?view=overview](https://data.mds.gbrrestoration.org/dataset/102.100.100/653201?view=overview)

## Config File Path Descriptions

### Domain

Paths to different ADRIA domains or historic input files
- `reefmod_domain` : Domain found on teams called `limited_reefmod_domain` in ADRIA domain folder.
- `rme_domain` : Path to ReefModEngine

### Geospatial

Paths to geopackages, shapefiles and geospatial location data.
- `canonical_path` : Path to canonical geopackage
- `ltmp_shp` : Path to LTMP regional shape files
- `classification_path` : CSV file contains location classes in the same order as the ADRIA domain.
- `bioregion_group_gpkg` : <!-- TODO -->

### Observations

Paths to data used as the target observation data for calibration.
- `manta_tow_path` : NetCDF containing target mean and standard deviation for location classes not individual
locations. Contained in `ltmp_data` directory.
- `ltmp_reef_data` : Geopackage containing target data for individual locations. Contained in `ltmp_data` directory.
- `composition_netcdf` : NetCDF containing coral composition for each ADRIA functional group at each ltmp
photogrammetry location. Contained in the `ltmp_data` directory.
- `ltmp_modelled_obs` : Path to modelled regional coral cover data based on LTMP
  observations

### Initialisation
- `init_cover_filepath` : Data containing calibrated initial cover. Must be loaded into domain as follows.
- `init_guess_filepath` : Optional file name for initial guess.

Initial cover must be loaded as follows.
```julia
init_cover = deserialize(init_cover_fn)
construct_cover!(dom, init_cover, location_classification.consecutive_classification)
```

### Outputs
- `out_dir` : Directory to save results plots, intermediate progress reports and final
  calibration results
  `result_filename` : Name of file to save calibrated results to, relative to `out_dir`

## Running scripts

The scripts in the `src` directory should be executed in the order according to the file
name numbering (excluding `2_calibration.jl`).

#### Running calibration and analysis

```julia-repl
julia> include("1_setup.jl")
julia> include("2b_location_calibration.jl")
julia> include("3_regional_analysis.jl")
julia> include("4_location_analysis.jl")
```

#### Results analysis only.

The calibration parameters used for results analysis are read from the `result_filename`
file in the output directory, `out_dir`, that were defined in the calibration configuration file.
Any plots already existing in this directory that were created using the same scripts will
be **overwritten**.

```julia-repl
julia> include("1_setup.jl")
julia> include("1b_data_split.jl")
julia> include("3_regional_analysis.jl")
julia> include("4_location_analysis.jl")
```

## Changing Calibration Parameters

In order to reduce the number of files to be edited when making changes to which parameters
are calibrated and what bounds they have, all parameter bounds setup is confined to the
`common/param_bounds.jl` file.

The `setup_run` function accepts a vector of parameter values and constructs the domain,
scenarios and scale factors required to perform a model run.

## Useful Variables

First execute `1_setup.jl`

- `dom` : ADRIA reefmod domain.
- `rm_ltmp_taxa` : Observed taxa composition at target locations
- `raw_ltmp_reef_data` : Observed coral cover levels at target locations
- `TARGET_DOM_IDXS` : ADRIA domain row index for each target location
- `NORTH MASK`, `CENTRAL MASK` and `SOUTH MASK` : LTMP Region mask for ADRIA domains.
- `ltmp_north`, `ltmp_central` and `ltmp_south` : Regional LTMP data
- `ALL_LTMP_REEF` : All ltmp reef manta cover data
- `ALL_LTMP_REEF_IDXS` : ADRIA domain row index for each ltmp location.
- `location_classification` : CSV containing location classification of GBR-wide locations

## Calibration-specific branches

Both ADRIA and CoralBlox have been modified to facilitate the calibration process.

The `ADRIA.jl#takuya/calib` contains a run model function with an altered call signature:

`rs_raw = ADRIA.run_model(dom, scens, scale_factors, target_loc_idxs, bleaching_threshold)`

where `scale_factors` is a matrix of shape `[taxa x n_target_locs x n_parameter_scaled]`.
The scale factors apply a different scaling coefficient for each location.

The bleaching threshold parameter specifies the lower bound of the truncated normal
distribution used to calculate bleaching mortality.

The latest `CoralBlox.jl#main` contains an implementation of `linear_extension_scale_factors`
which accepts parameters for a single location. It allows parameter values for a specific
location to be changed:

`linear_extension_scale_factors(::Matrix{Float64}, ::Float64, ::Matrix{Float64}, ::Matrix{Float64}, ::Float64)`

In addition to the usual:

`linear_extension_scale_factors(::AbstractArray{Float64, 3}, ::AbstractVector{Float64}, ::AbstractMatrix{Float64}, ::AbstractMatrix{Float64}, ::AbstractVector{Float64})`

## Error Functions

### Coral Cover Error

`reef_error(cover)`

For each target location calculate the mean absolute exponential error. Then calculate the
average over the given locations.

Expects total cover of shape `[timesteps x location]`

### Coral Composition Error

`reef_taxa_error(cover)`

For each target location calculate the correlation between ADRIA-mod coral composition and
the target composition contained in the coral_composition netcdf.

Expects cover of shape `[timesteps x taxa x size_class x location]`

### Location Class Error

`class_error(cover)` ~ calculates the average number of standard deviations a class is away
from the target cover for a given location class. The 'standard deviations' is based on the
spread of ltmp locations for a given year with a set of locations (class).

Expects total cover of shape `[timesteps x location]`

### Result Analysis

#### results_analysis.jl
