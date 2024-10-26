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

```bash
# TODO: Detail which entries require relative vs absolute paths
[data_paths]
reefmod_domain = "<path to ReefMod dataset>"  # Required
rme_domain = "<path to RME dataset>"  # Required
ltmp_shp = "spatial_data/gbr_3Zone 2.shp"  # Required
ltmp_modelled_obs = "ltmp_data/modelled_brms.beta.ry.disp.csv"  # Required
canonical_path = "<path to canonical gpkg>"  # Required
classification_path = "spatial_data/location_classification_MPA.csv"  # location classes
manta_tow_path = "ltmp_data/manta_tow_mean_std.nc"  # target data for location classes
ltmp_reef_data = "ltmp_data/manta_tow_data_reef_lvl.gpkg"  # target data for ltmp locs
composition_netcdf = "ltmp_data/coral_composition.nc"
init_cover_fn = "spatial/init_cover.dat"
init_guess_fn = "coral_p_calib_last.dat"  # optional, just the filename
out_dir = "../outputs"  # path to output directory for intermediate calibration results
```

## Config File Path Descriptions

- `reefmod_domain` : Domain found on teams called `limited_reefmod_domain` in ADRIA domain folder.

- `classification_path` : CSV file contains location classes in the same order as the ADRIA domain.

- `manta_tow_path` : NetCDF containing target mean and standard deviation for location classes not individual
locations. Contained in `ltmp_data` directory.

- `ltmp_reef_data` : Geopackage containing target data for individual locations. Contained in `ltmp_data` directory.

- `composition_netcdf` : NetCDF containing coral composition for each ADRIA functional group at each ltmp
photogrammetry location. Contained in the `ltmp_data` directory.

- `init_guess_fn` : Optional file name for inital guess.

- `init_cover_fn` : Data containing calibrated initial cover. Must be loaded into domain as follows.

```julia
init_cover = deserialize(init_cover_fn)
construct_cover!(dom, init_cover, location_classification.consecutive_classification)
```

## Useful Variables

First execute `1_setup.jl`

 - `dom` : ADRIA reefmod domain.

 - `rm_ltmp_taxa` : Target taxa composition at target locations

 - `raw_ltmp_reef_data` : Target coral cover levels at target locations

 - `target_dom_idxs` : ADRIA domain row index for each target location

 - `NORTH MASK`, `CENTRAL MASK` and `SOUTH MASK` : LTMP Region mask for ADRIA domains.

 - `ltmp_north`, `ltmp_central` and `ltmp_south` : Regional LTMP data

 - `all_ltmp_reef` : All ltmp reef manta cover data

 - `all_ltmp_reef_idxs` : ADRIA domain row index for each ltmp location.

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
