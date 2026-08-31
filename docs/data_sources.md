# Data Sources

Every dataset the calibration and analysis pipelines read, with its origin, licence, and
the config key or code path that consumes it.

- **Bundled** datasets ship in this repository under `datasets/` and `src/**/data/`.
- **External** datasets must be supplied by the user via `config.toml`.

Derived AIMS Long-Term Monitoring Program (LTMP) products are prepared by
[Data-Gen-Calibration](https://github.com/DanTanAtAims/Data-Gen-Calibration) and the
[ReefMonitoring-DataClient](https://github.com/open-AIMS/ReefMonitoring-DataClient) web
client.

## External inputs

| Config key | Dataset | Source | Used for |
|---|---|---|---|
| `rme_domain` (required) | ReefMod Engine data package (ReefMod Engine v1.0.28; e.g. `rme_ml_2024_01_08`) | ReefMod-GBR model (Bozec et al. 2025), <https://github.com/ymbozec/REEFMOD.7.2_GBR>. Available on request (y.bozec@uq.edu.au) or via the RRAP M&DS Data Store | Domain geometry, reef list, spatial layers; base for the ADRIA RME domain |
| `canonical_path` (required) | Canonical GBR reef geopackage (`rrap_canonical_<timestamp>.gpkg`) | Built by [gbrrestoration/canonical-reefs](https://github.com/gbrrestoration/canonical-reefs). A pre-built geopackage is available from that repository | Overwrites reef classification/identity columns in the RME domain on load; reef↔GBRMPA-ID mapping in the disturbance-mortality pipeline |
| `out_dir` (required) | — | user-chosen | Calibration outputs (`results.dat`, `params/*.nc`, plots) |
| — | `CoralSea_GBR_coraltempv3p1_dhw_1985-2024-reefs.mat` | AIMS M&DS IS Store | Input to `src/historical_dhw_data_gen/`; only needed to regenerate `historical_dhw.nc` (gitignored) |

**Canonical geopackage licence.** The `canonical-reefs` code is MIT-licensed. The output
geopackage has no single stated licence; it inherits terms from its source layers (GBRMPA
features, GBRMPA bioregions and zones, the ReefMod v1.0.28 reef list, AIMS LTMP cover, 10 m
bathymetry, EcoRRAP sites, and others). One source layer — the COTS control data — is
CC-BY-NC 4.0 (non-commercial).

## Bundled — third party

### GBR Features (GBRMPA)

`datasets/spatial_data/Great_Barrier_Reef_Features.geojson`

- **Contents**: reef boundaries, QLD mainland, islands, cays, rocks and dry reefs.
- **Custodian**: Great Barrier Reef Marine Park Authority (GBRMPA).
- **Source**: eAtlas, metadata record `ac8e8e4f-fc0e-4a01-9c3d-f27e4a8fac3c`
  (<https://catalogue.eatlas.org.au/geonetwork/srv/api/records/ac8e8e4f-fc0e-4a01-9c3d-f27e4a8fac3c>);
  also on GBRMPA Reef GeoHub.
- **Licence**: CC-BY 4.0. Attribution: *"© Great Barrier Reef Marine Park Authority 2014.
  Updated data available at http://www.gbrmpa.gov.au/geoportal/."*
- **Used by**: `src/viz/viz.jl` (`GBRMPA_MAINLAND_POLYS`, map basemap) and the reef-to-coast
  haversine distance in the spatial-correlation analysis.

### LTMP regions

`datasets/spatial_data/gbr_3Zone 2.{shp,shx,dbf,prj}` (and an identical copy at
`initial_cover/data/gbr_3Zone 2.*`)

- **Contents**: the LTMP Northern / Central / Southern GBR reporting regions — 3 polygons
  spanning the GBR Marine Park extent, split into latitudinal bands (Northern|Central
  boundary at −15.40°S; the Central and Southern polygons overlap between roughly −19.9°S
  and −20.7°S). The attribute table holds only `FID`, so region identity is by row order
  (Northern, Central, Southern). Consuming code depends on this order.
- **Source**: an internal AIMS LTMP regional-boundary layer, not published in a catalogue.
  Contact the AIMS monitoring team (<monitoring@aims.gov.au>) for the authoritative version
  and the exact region definitions. `.prj` is WGS84 with the datum name removed.
- **Licence**: AIMS data, CC-BY 4.0 (<https://www.aims.gov.au/cc-attribution>).
- **Used by**: `src/viz/regional_data_loader.jl` (`ltmp_shp` config default, regional
  analysis plots); `initial_cover/1_setup.jl` (ecological-class construction).

### AIMS LTMP observations

Derived from AIMS LTMP data, served through the Reef Monitoring dashboard
(<https://apps.aims.gov.au/reef-monitoring/>; metadata record
`5bb9a340-4ade-11dc-8f56-00008a07204e`).

- **Licence**: CC-BY 4.0 (<https://www.aims.gov.au/cc-attribution>). These products are
  derivative; attribute as *"Based on data supplied by: Australian Institute of Marine
  Science."*

| File | Contents | Config key | Used by |
|---|---|---|---|
| `datasets/ltmp_data/manta_tow_data_reef_lvl.gpkg` | Reef-level manta-tow total coral cover, keyed to `RME_UNIQUE_ID` | `ltmp_reef_data` | Calibration/test target for total cover |
| `datasets/ltmp_data/coral_composition.nc` | Community composition per ADRIA functional group at photogrammetry sites | `composition_netcdf` | Calibration/test target for composition |
| `datasets/ltmp_data/manta_tow_mean_std.nc` | Manta-tow cover mean and std aggregated by location class | `manta_tow_path` | `initial_cover/` pipeline — location-class error term |
| `datasets/ltmp_data/disturbances.nc` | Disturbance events, presence/absence (see below) | — (hard-coded path) | `scripts/plot/viz_results.jl` diagnostic plots |
| `datasets/ltmp_data/modelled_brms.beta.ry.disp.csv` | Modelled regional coral cover (Northern/Central/Southern, with CIs); Logan & Emslie BRMS model | `ltmp_modelled_obs` | Regional analysis plots; initial-cover setup |

**`datasets/ltmp_data/disturbances.nc`** — one Int16 variable `layer` over `locations`
(204 reefs, 11-digit GBRMPA reef IDs) × `timesteps` (2008–2022) × `disturbances`
(`Unknown`, `Storm`, `Multiple`, `Bleaching`, `Disease`, `Cots`); `1` marks a disturbance of
that type recorded at that reef in that year. Used to overlay observed disturbances on the
per-reef diagnostic plots. Built from `ReefMonitoring.get_disturbances()` output; the
original build script is not in the repo. To regenerate: for each reef from
`ReefMonitoring.get_reef_info()`, take the `MANTA` sample-type rows of
`ReefMonitoring.get_disturbances(reef_name)`, set `layer = 1` per `(reef, year, type)` for
2008–2022 (API codes `u/s/m/b/d/c` → the six labels above), and save as a YAXArray dataset.

### NOAA thermal stress

`src/historical_dhw_data_gen/data/historical_dhw.nc`

- **Contents**: annual maximum Degree Heating Weeks per reef, 1985–2024 (each year covers
  1 July of the previous year to 30 June).
- **Source**: derived from the NOAA Coral Reef Watch CoralTemp v3.1 daily product
  (<https://coralreefwatch.noaa.gov/product/5km/>).
- **Licence**: NOAA CoralTemp — no restrictions; attribution requested.
- **Used by**: `historical_dhw_path`; historical DHW forcing attached to the domain. The
  calibration timeframe is subset from this at load time.

## Bundled — generated by this repository

| File | Produced by | Contents / licence |
|---|---|---|
| `src/historical_dhw_data_gen/data/historical_dhw.nc` | `src/historical_dhw_data_gen/` (from the NOAA `.mat`) | See NOAA thermal stress above |
| `.../historical_disturbance_mortality_rates/historical_disturbance_mortality_rates.nc` | `src/historical_disturbance_mortality_data_gen/` (semi-automated, from LTMP disturbance and cover data) | Storm/COTS mortality per functional group, 2008–2022. CC-BY 4.0 (derived from AIMS LTMP data); see its co-located README |
| `datasets/spatial_data/location_classification_MPA.csv` (also `initial_cover/data/`) | Data-Gen-Calibration classification step | Ecological class per reef; `classification_path` default |
| `datasets/spatial_data/init_cover.dat` | A previous calibration run | Serialised initial-cover array; `init_cover_filepath` default |

## Attribution

Reproduce in `NOTICE` and in any publication using these data:

> Reef and coastline geometries are from the *Great Barrier Reef Features* dataset
> (© Great Barrier Reef Marine Park Authority 2014; CC-BY 4.0), obtained via the eAtlas
> data portal (record `ac8e8e4f-fc0e-4a01-9c3d-f27e4a8fac3c`).
>
> Coral cover, community composition and disturbance data are based on data supplied by the
> Australian Institute of Marine Science (Long-Term Monitoring Program), CC-BY 4.0
> (<https://www.aims.gov.au/cc-attribution>).
>
> Thermal-stress data are derived from the NOAA Coral Reef Watch CoralTemp v3.1 product.
>
> Domain and environmental forcing data are from the ReefMod-GBR model (Bozec et al.).
