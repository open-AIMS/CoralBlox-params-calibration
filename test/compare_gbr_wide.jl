"""
The script should produce coral cover plots using the 'global-calib' ADRIA branch to compare
with the older 'takuya/calib' branch to make sure the implementation align.
"""

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

scens = ADRIA.param_table(dom)

grouped_bioregion_gpkg = GDF.read(
    "C:/Users/dtan/repos/Data-Gen-Calibration/Outputs/domain_bioregions.gpkg"
)
dom.loc_data[!, :GROUPED_BIOREGION] .= grouped_bioregion_gpkg.ASSIGNED_BIOREGION

unique_biogroups = unique(dom.loc_data.GROUPED_BIOREGION)

# names of dimension [taxa ⋅ param_type ⋅ biogroup]
scale_factor_names = ADRIA.generate_scale_factor_names(unique_biogroups)
new_col_names = ADRIA.scale_factor_array_to_vec(scale_factor_names)

insertcols!(scens, (new_col_names .=> Ref([1.0]))...)
scens[!, new_col_names] .= 1.0

# names of dimension [biogroups ⋅ params]
growth_accel_names = ADRIA.generate_growth_accel_names(unique_biogroups)
new_col_names = ADRIA.accel_params_array_to_vec(growth_accel_names)

insertcols!(scens, (new_col_names .=> Ref([1.0]))...)
growth_accel_params = ones(Float64, length(unique_biogroups), 3)
growth_accel_params[:, 2] .= 0.0
scens[!, new_col_names]   .= ADRIA.accel_params_array_to_vec(growth_accel_params)'

rs_raw = ADRIA.run_model(dom, scens[1, :])

using CairoMakie

k_area = ADRIA.site_k_area(dom)

relative_cover = dropdims(sum(rs_raw.raw .* reshape(k_area, (1, 1, 3806)), dims=(2, 3)) ./ sum(k_area), dims=(2, 3))

f = Figure()
Axis(f[1, 1], title="test scenario run (ADRIA/global-calib)", xlabel="Year", ylabel="Relative Cover")
lines!(2008:2022, relative_cover)
save("../test/relative_cover_global_calib.png", f)

"""
Mimic the target locations of the older branch 'takuya/calib' by assigning each target
location their own bioregion and all other locations the same bioregion. Set the scale
factors of all other locations to be one and the growth acceleration parameters to be the
mean of the target location parameters. This should produce equivalent outputs.
"""
target_loc_idxs = [54, 876, 234, 1234]

dom.loc_data[:, :GROUPED_BIOREGION] .= 1
for (idx, loc_idx) in enumerate(target_loc_idxs)
    dom.loc_data[loc_idx, :GROUPED_BIOREGION] = idx + 1
end

unique_biogroups = sort(unique(dom.loc_data.GROUPED_BIOREGION))

# names of dimension [taxa ⋅ param_type ⋅ biogroup]
scale_factor_names = ADRIA.generate_scale_factor_names(unique_biogroups)
scale_factors = ones(Float64, size(scale_factor_names)...)

# The following scale factors are not ecologically meaningful and are for testing only.
# linear extension scale factors
scale_factors[:, LIN_EXT_SCALE_IDX, 2] .= 1.2
scale_factors[:, LIN_EXT_SCALE_IDX, 3] .= 1.4
scale_factors[:, LIN_EXT_SCALE_IDX, 4] .= 0.8
scale_factors[:, LIN_EXT_SCALE_IDX, 5] .= 0.6

# Background mortality scale factors
scale_factors[:, MORTALITY_SCALE_IDX, 2] .= 0.95
scale_factors[:, MORTALITY_SCALE_IDX, 3] .= 0.9
scale_factors[:, MORTALITY_SCALE_IDX, 4] .= 0.85
scale_factors[:, MORTALITY_SCALE_IDX, 5] .= 0.8

# Reshape parameter names and parameters to be used for assignment in scenario dataframe
new_col_names = ADRIA.scale_factor_array_to_vec(scale_factor_names)
scale_factors = ADRIA.scale_factor_array_to_vec(scale_factors)

scens = ADRIA.param_table(dom)
insertcols!(scens, (new_col_names .=> Ref([1.0]))...)
scens[!, new_col_names] .= scale_factors'

# The following growth accerlation parameters are not ecologically meaningful and are for
# testing purposes only.

# names of dimension [biogroups ⋅ params]
growth_accel_names = ADRIA.generate_growth_accel_names(unique_biogroups)
growth_accel_params = ones(Float64, size(growth_accel_names)...)

# Define the parametrisation of growth acceleration for the non-target locations
growth_accel_params[1, STEEPNESS_IDX] = -20.0
growth_accel_params[1, MIDPOINT_IDX]  = 0.3
growth_accel_params[1, HEIGHT_IDX]    = 0.4

# Define the parametrisation of growth acceleration for the target locations
growth_accel_params[2, STEEPNESS_IDX] = -25.0
growth_accel_params[3, STEEPNESS_IDX] = -30.0
growth_accel_params[4, STEEPNESS_IDX] = -15.0
growth_accel_params[5, STEEPNESS_IDX] = -10.0

growth_accel_params[2, HEIGHT_IDX] = 0.5
growth_accel_params[3, HEIGHT_IDX] = 0.6
growth_accel_params[4, HEIGHT_IDX] = 0.3
growth_accel_params[5, HEIGHT_IDX] = 0.2

growth_accel_params[2, MIDPOINT_IDX] = 0.35
growth_accel_params[3, MIDPOINT_IDX] = 0.4
growth_accel_params[4, MIDPOINT_IDX] = 0.25
growth_accel_params[5, MIDPOINT_IDX] = 0.2

@assert (
    mean(growth_accel_params[2:5, STEEPNESS_IDX]) == growth_accel_params[1, STEEPNESS_IDX]
) "Mean of target location steepness param not equal to non-target location param"
@assert (
    mean(growth_accel_params[2:5, HEIGHT_IDX]) == growth_accel_params[1, HEIGHT_IDX]
) "Mean of target location height param not equal to non-target location param"
@assert (
    mean(growth_accel_params[2:5, MIDPOINT_IDX]) == growth_accel_params[1, MIDPOINT_IDX]
) "Mean of target location midpoint param not equal to non-target location param"

new_col_names = ADRIA.accel_params_array_to_vec(growth_accel_names)
growth_accel_params = ADRIA.accel_params_array_to_vec(growth_accel_params)

insertcols!(scens, (new_col_names .=> Ref([1.0]))...)
scens[!, new_col_names] .= growth_accel_params'

rs_raw = ADRIA.run_model(dom, scens[1, :])

using CairoMakie

k_area = ADRIA.site_k_area(dom)

relative_cover = dropdims(sum(rs_raw.raw .* reshape(k_area, (1, 1, 3806)), dims=(2, 3)) ./ sum(k_area), dims=(2, 3))

f = Figure()
Axis(f[1, 1], title="test scenario run with model changes (ADRIA/global-calib)", xlabel="Year", ylabel="Relative Cover")
lines!(2008:2022, relative_cover)
save("../test/relative_cover_global_calib_model_changes.png", f)

loc_save_dir = joinpath("..", "test", "locations")
mkpath(loc_save_dir)

for loc_idx in target_loc_idxs
    local coral_cover = dropdims(sum(rs_raw.raw[:, :, loc_idx], dims=2), dims=2)
    local f = Figure()
    Axis(f[1, 1], title="location $(loc_idx) with model changes (ADRIA/global-calib)", xlabel="Year", ylabel="Relative Cover")
    lines!(2008:2022, coral_cover)
    save(joinpath(loc_save_dir, "new_branch_location_$(loc_idx).png"), f)
end
