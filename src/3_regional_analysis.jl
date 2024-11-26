using ADRIA: RESULTS
include("./common/common.jl")
include("./common/cover_construction.jl")
include("./common/plotting/plotting.jl")
include("1_setup.jl")
include("common/param_bounds.jl")

init_cover_fn = "C:/Users/dtan/data/init_cover.dat"
init_cover = deserialize(init_cover_fn)

construct_cover!(
    dom, init_cover, location_classification.consecutive_classification
)

# ----- LOAD CALIBRATED RESULTS -----

coral_param_fn = joinpath(OUT_DIR, RESULT_FN)
coral_params = deserialize(coral_param_fn)

# Load values into scenario dataframe

dom, scen, growth_acc_params, scale_factors = setup_run(
    dom,
    coral_params;
    param_names=CORAL_PARAM_NAMES,
    param_idxs=PARAM_IDXS,
    loc_idxs=TARGET_DOM_IDXS
)

rs_raw = ADRIA.run_model(dom, scen[1, :], scale_factors, growth_acc_params, TARGET_DOM_IDXS)

s_rac = (dropdims(sum(rs_raw.raw, dims=2), dims=2) .* site_k_area(dom)') ./ loc_area(dom)'

ref_years = START_YEAR:END_YEAR

north_res = ADRIA.DataCube(s_rac[:, NORTH_MASK, :]; timesteps=ref_years, sites=1:count(NORTH_MASK), scenarios=1:1)
central_res = ADRIA.DataCube(s_rac[:, CENTRAL_MASK, :]; timesteps=ref_years, sites=1:count(CENTRAL_MASK), scenarios=1:1)
south_res = ADRIA.DataCube(s_rac[:, SOUTH_MASK, :]; timesteps=ref_years, sites=1:count(SOUTH_MASK), scenarios=1:1)

f = plot_all_regions(
    north_res, ltmp_north, central_res, ltmp_central, south_res, ltmp_south
)

save_dir = OUT_DIR

mkpath(save_dir)

save(joinpath(save_dir, "locs_reg.png"), f)

f = taxa_cover_proportions(rs_raw.raw)

save(joinpath(save_dir, "locs_taxa_cov.png"), f)

f = taxa_population_proportions(rs_raw.raw)

save(joinpath(save_dir, "locs_taxa_pop.png"), f)

f = temporal_size_class_proportions(rs_raw.raw)

save(joinpath(save_dir, "locs_size.png"), f)

mkpath("$(save_dir)/loc_plots")

for i in 1:296
    location_comparison(rs_raw.raw, i, "$(save_dir)/loc_plots")
end
