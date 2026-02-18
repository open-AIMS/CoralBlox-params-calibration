# ADRIA-CoralBlox calibration

This repo is used to calibrate some parameters of ADRIA/CoralBlox model to match, as much as
possible, historic LTMP (Long-Term Monitoring Program) data for the GBR (Great Barrier Reef)
between years 2008 and 2022. The parameters calibrated here were:

- `mb_rate`: Base mortality rate for each functional group and size class
- `linear_extension`: Base linear extension for each functional group and size class
- `biogroup_linear_extension`: Scale factor to be applied to linear extensions of all
functional groups and size classes of a spatial group. We are now using `spatial_group`
instead of `biogroup` and need to update that in this project.
- `biogroup_mb_rate`: Scale factor to be applied to base mortality rates of all
functional groups and size classes of a spatial group. We are now using `spatial_group`
instead of `biogroup` and need to update that in this project.
- `growth_accel_steepness`, `growth_accel_height` and `growth_accel_midpoint`: Used as
parameters of a growth acceleration function.

The equation used for the growth acceleration is:

```
height / (1 + exp(-steepness * (available_space - midpoint))) + 1.0
```

As part of this calibration, two datasets were generated, one with historical DHW values
and another with historical cyclone/storm/COTS mortality rates. Further details on how these
 can be found in sections below.

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
rme_domain = <path to RME dataset>

[Geospatial]
canonical_path = <path to canonical gpkg>
ltmp_shp = "../datasets/spatial_data/gbr_3Zone 2.shp"                             # Optional
classification_path = "../datasets/spatial_data/location_classification_MPA.csv"  # Optional

[Observations]
manta_tow_path = "../datasets/ltmp_data/manta_tow_mean_std.nc"                    # Optional
ltmp_reef_data = "../datasets/ltmp_data/manta_tow_data_reef_lvl.gpkg"             # Optional
composition_netcdf = "../datasets/ltmp_data/coral_composition.nc"                 # Optional
ltmp_modelled_obs = "../datasets/ltmp_data/modelled_brms.beta.ry.disp.csv"        # Optional

[Initialisation]
init_cover_filepath = "../datasets/spatial_data/init_cover.dat"                   # Optional
init_guess_filepath = ""                                                          # Optional
ecorrap_param_filepath = <path to interped_vals.nc>

[Outputs]
out_dir = <path to output dir>
result_filename = <path to results.dat>                               # relative to out_dir
```

## Datasets

### Domains

Paths to different ADRIA domains or historic input files
- `rme_domain` : Path to ReefModEngine dir "rme_ml_2025_06_05"

### Geospatial

Paths to geopackages, shapefiles and geospatial location data.
- `canonical_path` : Path to canonical geopackage
- `ltmp_shp` : Path to LTMP regional shape files
- `classification_path` : CSV file contains location classes in the same order as the ADRIA
domain.

### Observations

Paths to data used as the target observation data for calibration.
- `manta_tow_path` : NetCDF containing target mean and standard deviation for location
classes not individual locations. Contained in `ltmp_data` directory.
- `ltmp_reef_data` : Geopackage containing target data for individual locations. Contained
in `ltmp_data` directory.
- `composition_netcdf` : NetCDF containing coral composition for each ADRIA functional group
at each ltmp photogrammetry location. Contained in the `ltmp_data` directory.
- `ltmp_modelled_obs` : Path to modelled regional coral cover data based on LTMP
  observations. This was developed by Murray Logan and Mike Emslie.

### Initialisation

- `init_cover_filepath` : Data containing calibrated initial cover. Must be loaded into
domain as follows.
- `init_guess_filepath` : Optional file name for initial guess.
- `ecorrap_param_filepath` : <!--TO DO-->

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
julia> include("01_a_setup.jl")
julia> include("01_b_data_split.jl")
julia> include("02_location_calibration.jl")
```

#### Results analysis only.

The calibration parameters used for results analysis are read from the `result_filename`
file in the output directory, `out_dir`, that were defined in the calibration configuration
file. Any plots already existing in this directory that were created using the same scripts
will be **overwritten**.

```julia-repl
julia> include("01_a_setup.jl")
julia> include("01_b_data_split.jl")
julia> include("03_result_analysis.jl")
julia> include("04_location_analysis.jl")
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
- `NORTH MASK`, `CENTRAL MASK` and `SOUTH MASK` : LTMP Region mask for ADRIA domains
- `ltmp_north`, `ltmp_central` and `ltmp_south` : Regional LTMP data
- `ALL_LTMP_REEF` : All ltmp reef manta cover data
- `ALL_LTMP_REEF_IDXS` : ADRIA domain row index for each ltmp location
- `location_classification` : CSV containing location classification of GBR-wide locations

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

## Historical DHW

The script `src/historical_dhw_data_gen/01_generate_historical_dhw_data.jl` takes a file
containing observed NOAA maximum DHWs at each of the reefs, selects the timeframe used in
this calibration (2008-2022) and saves the result as a NetCDF file in
`src/historical_dhw_data_gen/data/historical_dhw.nc`. Before running this script, the input
file `CoralSea_GBR_coraltempv3p1_dhw_1985-2024-reefs.mat` must be placed at
`src/historical_dhw_data_gen/data/`. This file can be found in the M&DS IS Store.

The data in the input file is stored in a MATLAB structure array ‘R’. From this file, we
have extracted the variables:

- `years_dhw` : 1985 to 2024, where 1985 only contains from early 1985 as the SST dataset
only starts from 1 Jan. Otherwise each year is for the maximum DHW from 1 July year-1 to 30
June year, with the final year being from 1 July 2023 to 30 June 2024.
- `dhw_max` : observed maximum DHWs at each reef

## Historical cyclone/storms/COTS mortality rates

Scripts in `src/historical_disturbance_mortality_data_gen` are used to generate the
`historical_disturbance_mortality_rates.nc` dataset. This process is semi-automated. After
running script 01, the csv files created (`template_disturbance_years_manta.csv` and `template_disturbance_years_transect.csv`) should be manually filled with the years before
and after each disturbance with data from the Reef Monitoring Dashboard
(https://apps.aims.gov.au/reef-monitoring/reefs). These csv files also need to be renamed to ("disturbance_years_manta.csv" and "disturbance_years_transect.csv") before the next script
is run. To our knowledge, there is no automatic way of determining which year is before and
after a disturbance.

More details about this dataset can be found in `src\historical_disturbance_mortality_data_gen\data\historical_disturbance_mortality_rates\README.md`