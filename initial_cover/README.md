# ADRIA-CoralBlox Initial Cover calibration

Code taken from an [older version](https://github.com/ConnectedSystems/ltmp_calibration/tree/92894d434d95dc4b9b73ad60c8ef8a4e39bcba7c) of this repository

This repository contains the initial attempts at calibrating the starting
coral cover of ADRIA-CoralBlox.

The classification subdirectory of [Data-Gen-Calibration](https://github.com/DanTanAtAims/Data-Gen-Calibration/tree/main/src/classification)
classifies locations in the canonical geopackage by aggregating depth, turbidty, and wave
activity into thirds and then also groups locations by shelf position. This calibration
process seeks to approximate coral cover levels, coral composition and size class
distribution for each location classification. This data is contained in the
`location_classication_MPA.csv` file.

Size class distribution is assumed to be exponentially distributed and the calibration
process calibrates the rate parameter for each size class.

The final calibrated params will be saved in `Outputs/coral_cover.dat`

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
ltmp_shp = "data/gbr_3Zone 2.shp"
canonical_path = "<path to canonical geopackage>"
classification_path = "data/location_classification_MPA.csv"
```

**The shape files and classification files can be found in the data subdirectory.**
