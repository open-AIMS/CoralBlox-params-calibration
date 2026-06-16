using GeoDataFrames
using MAT
using NetCDF
using YAXArrays
using TOML

isdefined(Main, :CONFIG) || (global CONFIG = TOML.parsefile("config.toml"))
isdefined(Main, :GEOSPATIAL_CONFIG) || (global GEOSPATIAL_CONFIG = CONFIG["calibration"]["geospatial"])
isdefined(Main, :CANONICAL_PATH) || (global CANONICAL_PATH = GEOSPATIAL_CONFIG["canonical_path"])

begin
    canonical_gpkg = GeoDataFrames.read(CANONICAL_PATH)

    dhw_mat_path = "data/CoralSea_GBR_coraltempv3p1_dhw_1985-2024-reefs.mat"

    vars = matread(dhw_mat_path)

    # Relevant variables holding the years and DHW for each reef on Canonical Reefs
    data = vars["R"]["dhw_max"]'

    n_timesteps, n_locations = size(data)

    timesteps = Int64.(dropdims(vars["R"]["years_dhw"], dims=2))
    locations = canonical_gpkg.RME_UNIQUE_ID

    ax = (
        Dim{:timesteps}(first(timesteps):last(timesteps)),
        Dim{:location}(locations),
        Dim{:scenarios}(1:2)
    )

    dhw_scens = YAXArray(ax, zeros(n_timesteps, n_locations, 2))
    dhw_scens .= data

    # Save as NetCDF
    output_path = "data/dhw_scens.nc"
    savedataset(Dataset(; dhw_scens), path=output_path, driver=:netcdf)
end