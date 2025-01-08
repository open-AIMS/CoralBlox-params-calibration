using Infiltrator, Revise
using ADRIA

using CSV, DataFrames
import GeoDataFrames as GDF

include("../src/1_setup.jl")

LIN_EXT_SCALE_IDX = 1
MORTALITY_SCALE_IDX = 2

STEEPNESS_IDX = 1
HEIGHT_IDX = 2
MIDPOINT_IDX = 3

coefs = ones(Float64, 5, 1, 3)
growth_accel_params = ones(Float64, 3, 1)
growth_accel_params[2, 1] = 0.5
cloc_idxs = [1]

scens = ADRIA.param_table(dom)
rs_raw = ADRIA.run_model(dom, scens[1, :], coefs, growth_accel_params, cloc_idxs)

using CairoMakie

k_area = ADRIA.site_k_area(dom)

relative_cover = dropdims(sum(rs_raw.raw .* reshape(k_area, (1, 1, 3806)), dims=(2, 3)) ./ sum(k_area), dims=(2, 3))

f = Figure()
Axis(f[1, 1], title="test scenario run (ADRIA/takuya/calib)", xlabel="Year", ylabel="Relative Cover")
lines!(2008:2022, relative_cover)
save("../test/relative_cover_old.png", f)

cloc_idxs = [54, 876, 234, 1234]

# Coefficients should match coefficients found in 'test_gbr_wide.jl'
coefs = ones(Float64, 5, 4, 3)

coefs[:, 1, LIN_EXT_SCALE_IDX] .= 1.2
coefs[:, 2, LIN_EXT_SCALE_IDX] .= 1.4
coefs[:, 3, LIN_EXT_SCALE_IDX] .= 0.8
coefs[:, 4, LIN_EXT_SCALE_IDX] .= 0.6

coefs[:, 1, MORTALITY_SCALE_IDX] .= 0.95
coefs[:, 2, MORTALITY_SCALE_IDX] .= 0.9
coefs[:, 3, MORTALITY_SCALE_IDX] .= 0.85
coefs[:, 4, MORTALITY_SCALE_IDX] .= 0.8

growth_accel_params = ones(Float64, 3, 4)

growth_accel_params[STEEPNESS_IDX, 1] = -25.0
growth_accel_params[STEEPNESS_IDX, 2] = -30.0
growth_accel_params[STEEPNESS_IDX, 3] = -15.0
growth_accel_params[STEEPNESS_IDX, 4] = -10.0

growth_accel_params[HEIGHT_IDX, 1] = 0.5
growth_accel_params[HEIGHT_IDX, 2] = 0.6
growth_accel_params[HEIGHT_IDX, 3] = 0.3
growth_accel_params[HEIGHT_IDX, 4] = 0.2

growth_accel_params[MIDPOINT_IDX, 1] = 0.35
growth_accel_params[MIDPOINT_IDX, 2] = 0.4
growth_accel_params[MIDPOINT_IDX, 3] = 0.25
growth_accel_params[MIDPOINT_IDX, 4] = 0.2

scens = ADRIA.param_table(dom)
rs_raw = ADRIA.run_model(dom, scens[1, :], coefs, growth_accel_params, cloc_idxs)

k_area = ADRIA.site_k_area(dom)

relative_cover = dropdims(sum(rs_raw.raw .* reshape(k_area, (1, 1, 3806)), dims=(2, 3)) ./ sum(k_area), dims=(2, 3))

f = Figure()
Axis(f[1, 1], title="test scenario run with model changes (ADRIA/takuya/calib)", xlabel="Year", ylabel="Relative Cover")
lines!(2008:2022, relative_cover)
save("../test/relative_cover_old_model_changes.png", f)

loc_save_dir = joinpath("..", "test", "locations")
mkpath(loc_save_dir)

for loc_idx in cloc_idxs
    local coral_cover = dropdims(sum(rs_raw.raw[:, :, loc_idx], dims=2), dims=2)
    local f = Figure()
    Axis(f[1, 1], title="location $(loc_idx) with model changes (ADRIA/takuya/calib)", xlabel="Year", ylabel="Relative Cover")
    lines!(2008:2022, coral_cover)
    save(joinpath(loc_save_dir, "old_branch_location_$(loc_idx).png"), f)
end
