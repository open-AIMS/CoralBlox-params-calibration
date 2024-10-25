# ADRIA-CoralBlox calibration

## Setup

Instantiate environment to install required packages.

```julia
]instantiate
```

Add ADRIA using `dev`.

Here, it is assumed this repo lives inside the sandbox directory.
Otherwise, use the absolute path to the ADRIA repository.

```julia
] dev ../..
```

Create a `calib_config.toml` file with the following entries:

```bash
[data_paths]
reefmod_domain = "<path to ReefMod dataset>"
rme_domain = "<path to RME dataset>"
ltmp_shp = "spatial_data/gbr_3Zone 2.shp"
canonical_path = "<path to canonical gpkg>"
classification_path = "spatial_data/location_classification_MPA.csv" # location classes
manta_tow_path = "ltmp_data\\manta_tow_mean_std.nc" # target data for location classes
ltmp_reef_data = "ltmp_data\\manta_tow_data_reef_lvl.gpkg" # target data for ltmp locs
composition_netcdf = "ltmp_data\\coral_composition.nc"
init_cover_fn = "spatial\\init_cover.dat"
init_guess_fn = "<path to initial guess>" # optional
```

## Config File Path Descriptions

 - `reefmod_domain`

 Domain found on teams called `limited_reefmod_domain` in ADRIA domain folder.

 - `classification_path`

 CSV file contains location classes in the same order as the ADRIA domain.

 - `manta_tow_path`

 NetCDF containing target mean and standard deviation for location classes not individual
locations. Contained in `ltmp_data` directory.

 - `ltmp_reef_data`

 Gpkg containing target data for individual locations. Contained in `ltmp_data` directory.

 - `composition_netcdf`

 NetCDF containing coral composition for each ADRIA functional group at each ltmp
photogrammetry location. Contained in the `ltmp_data` directory.

 - `init_cover_fn`

 Data containing calibrated initial cover. Must be loaded into domain as follows.

```julia
init_cover = deserialize(init_cover_fn)
construct_cover!(dom, init_cover, location_classification.consecutive_classification)
```

 - `init_guess_fn`

 Optional file name for inital guess.

## ADRIA  Branch

Checkout `takuya/calib` for compatability.

The run model function has an altered call signature.

`rs_raw = ADRIA.run_model(dom, scens, scale_factors, target_loc_idxs, bleaching_threshold)`

where `scale_factors` is a matrix of shape `[taxa x n_target_locs x n_parameter_scaled]`.
The scale factors apply a different scaling coefficient for each location.

The bleaching threshold parameter specifies the lower bound of the truncated normal
disitrbution used to calculate bleaching mortality.

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


