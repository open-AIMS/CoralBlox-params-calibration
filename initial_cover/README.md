# ADRIA-CoralBlox Initial Cover Calibration

Code taken from an [older version](https://github.com/ConnectedSystems/ltmp_calibration/tree/92894d434d95dc4b9b73ad60c8ef8a4e39bcba7c) of this repository.

## Motivation

ADRIA requires an initial coral cover for every location in the domain (~3000 GBR reefs),
but historical observations from the Long-Term Monitoring Program (LTMP) exist for only a
subset of those reefs. This calibration process estimates initial cover for the unobserved
locations. At locations where historical data exists, that data is used directly (see below).

## Approach

### Location classification

The [classification subdirectory of Data-Gen-Calibration](https://github.com/DanTanAtAims/Data-Gen-Calibration/tree/main/src/classification)
assigns each reef in the canonical geopackage to an ecological class by aggregating depth,
turbidity, and wave energy into tertiles, then further stratifying by shelf position and LTMP
region. The resulting per-location class labels are stored in
`datasets/spatial_data/location_classification_MPA.csv` (column `consecutive_classification`).

All reefs in the same class are assumed to share the same initial cover profile. The
parameters of that profile — total cover, taxonomic composition weights, and a size-class
exponential rate per taxon — are calibrated jointly alongside the coral growth and mortality
parameters. Size-class distributions are modelled as squared-exponential distributions
parameterised by a rate λ per functional group.

### How observed and unobserved locations differ

Cover is assigned in two phases:

**Phase 1** stamps every location with the cover profile of its class. This is the only
cover assignment for locations with no historical data.

**Phase 2** overrides Phase 1 at locations where observations exist:

- *LTMP manta tow sites*: total cover is rescaled to match the earliest available
  observed cover; the composition shape from Phase 1 is preserved.
- *Photogrammetry sites*: both taxonomic composition and size-class distribution are
  replaced with observed values; only the total cover magnitude is retained from Phase 1.

The final calibrated cover profiles are saved to `Outputs/coral_cover.dat`.

## Setup

Instantiate the environment to install required packages.

```julia
]instantiate
```

Add ADRIA using `dev`. Here it is assumed this repo lives inside the sandbox directory;
otherwise use the absolute path to the ADRIA repository.

```julia
] dev ../..
```

Create a `calib_config.toml` file with the following entries:

```toml
[data_paths]
reefmod_domain = "<path to ReefMod dataset>"
rme_domain = "<path to RME dataset>"
ltmp_shp = "data/gbr_3Zone 2.shp"
canonical_path = "<path to canonical geopackage>"
classification_path = "data/location_classification_MPA.csv"
```

The shape files and classification file can be found in the `data` subdirectory.
