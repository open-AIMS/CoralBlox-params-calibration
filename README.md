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
canonical_path = "spatial_data/canonical_gbr_2024-04-23.gpkg"
```
