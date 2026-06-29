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

Copy `config.toml.template` to `config.toml` at the repo root and fill in your
machine-specific paths. The template documents all fields inline, including which are
required and what the relative-path defaults are for optional fields.

## Running scripts

The scripts in the `src` directory should be executed in the order according to the file
name numbering (excluding `2_calibration.jl`).

#### Running calibration and analysis

```julia-repl
julia> include("src/01_a_setup.jl")
julia> include("src/01_b_data_split.jl")
julia> include("src/02_location_calibration.jl")
```

#### Results analysis only.

The calibration parameters used for results analysis are read from the `result_filename`
file in the output directory, `out_dir`, that were defined in the calibration configuration
file. Any plots already existing in this directory that were created using the same scripts
will be **overwritten**.

```julia-repl
julia> include("src/01_a_setup.jl")
julia> include("src/01_b_data_split.jl")
julia> include("src/03_result_analysis.jl")
julia> include("src/04_location_analysis.jl")
```

## Initial Coral Cover

ADRIA requires an initial coral cover for every location in the domain. Because historical
coral cover observations (LTMP manta tow and photogrammetry data) are only available for a
subset of GBR reefs, a two-phase assignment is used:

### Phase 1 — Classification-based cover (all locations)

Every location is assigned an initial cover profile derived from its ecological class. Classes
are defined by aggregating depth, turbidity, wave energy, shelf position, and LTMP region into
a mixed-radix index (see `datasets/spatial_data/location_classification_MPA.csv`, produced by
[Data-Gen-Calibration](https://github.com/DanTanAtAims/Data-Gen-Calibration/tree/main/src/classification)).
All locations sharing the same class receive the same cover profile. The parameters of this
profile (total cover, taxonomic composition weights, and size-class exponential rate per
taxon) are jointly calibrated alongside the coral growth and mortality parameters.

### Phase 2 — Observation-based override (locations with historical data)

After Phase 1, cover is overwritten at locations where observations exist:

- *LTMP manta tow sites*: the total cover is rescaled to match the earliest non-missing
  observed cover, preserving the composition shape from Phase 1.
- *Photogrammetry sites*: both the taxonomic composition and size-class distribution are
  replaced with observed values; only the total cover magnitude is retained from Phase 1.

The index vectors `ltmp_cover_to_domain` and `composition_to_domain` inside each
`LocationDataStore` (`CALIBRATION_STORE`, `VALIDATION_STORE`, `COMBINED_STORE`) are the sole
record of which domain locations have observed data. Everything outside these index sets
retains its Phase 1 cover unchanged.

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

## Parameter Extraction

After calibration is complete, `src/common/params_extraction.jl` unpacks the flat
calibrated-parameter vector into labelled NetCDF files suitable for use in ADRIA.

**Prerequisites** — the following must be defined (run `01_a_setup.jl` first) and a
completed calibration result must exist at `joinpath(OUT_DIR, RESULT_FN)`:

- `dom` — ADRIA domain
- `OUT_DIR`, `RESULT_FN` — output directory and result filename (from `config.toml`)
- `INIT_COVER_PATH` — path to the serialised initial-cover sample

```julia-repl
julia> include("src/01_a_setup.jl")
julia> include("src/common/params_extraction.jl")
```

This writes two files to `{out_dir}/params/`:

- **`calibrated_params.nc`** — Calibrated coral parameters as a YAXArray `Dataset` with
  labelled `functional_group`, `size_class`, `cb_calib_group`, and `accel_param` axes.
  Variables: `linear_extension`, `mb_rate`, `dist_mean`, `linear_extension_scale`,
  `mb_rate_scale`, `growth_acceleration`.
- **`historic_init_cover.nc`** — GBR-wide historical initial coral cover with
  observation-based overrides applied at LTMP and photogrammetry sites
  (see [Initial Coral Cover](#initial-coral-cover)).

> **Note:** The `params/` subdirectory is created automatically if it does not exist.
> Existing files are overwritten without warning.

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
