using ADRIA: RESULTS
include("./common/common.jl")
include("./common/cover_construction.jl")
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
    loc_idxs=target_dom_idxs
)

rs_raw = ADRIA.run_model(dom, scens[1, :], scale_factors, growth_acc_params, target_dom_idxs)

s_rac = (dropdims(sum(rs_raw.raw, dims=2), dims=2) .* site_k_area(dom)') ./ loc_area(dom)'

ref_years = start_year:end_year
north_mean_cover = zeros(size(rs_raw.raw, 1))
center_mean_cover = zeros(size(rs_raw.raw, 1))
south_mean_cover = zeros(size(rs_raw.raw, 1))
comp_years_north = (ref_years .∈ [ltmp_north[ltmp_north.Year .>= start_year, :Year]])
comp_years_center = (ref_years .∈ [ltmp_central[ltmp_central.Year .>= start_year, :Year]])
comp_years_south = (ref_years .∈ [ltmp_south[ltmp_south.Year .>= start_year, :Year]])

for ts in axes(rs_raw.raw, 1)
    rac = vec(sum(rs_raw.raw[ts, :, :], dims=1) .* site_k_area(dom)') ./ loc_area(dom)
    north_mean_cover[ts] = mean(rac[NORTH_MASK])
    center_mean_cover[ts] = mean(rac[CENTRAL_MASK])
    south_mean_cover[ts] = mean(rac[SOUTH_MASK])
end

north_mean_cover = north_mean_cover[comp_years_north]
center_mean_cover = center_mean_cover[comp_years_center]
south_mean_cover = south_mean_cover[comp_years_south]

north_res = ADRIA.DataCube(s_rac[:, NORTH_MASK, :]; timesteps=ref_years, sites=1:count(NORTH_MASK), scenarios=1:1)
central_res = ADRIA.DataCube(s_rac[:, CENTRAL_MASK, :]; timesteps=ref_years, sites=1:count(CENTRAL_MASK), scenarios=1:1)
south_res = ADRIA.DataCube(s_rac[:, SOUTH_MASK, :]; timesteps=ref_years, sites=1:count(SOUTH_MASK), scenarios=1:1)


f = Figure(; Dict{Symbol,Any}(:size => (1600, 1600))...)
ax1 = plot_region(
    f,
    1,
    1,
    "North GBR",
    ltmp_north,
    north_res
)
ax2 = plot_region(
    f,
    1,
    2,
    "Central GBR",
    ltmp_central,
    central_res
)
ax3 = plot_region(
    f,
    2,
    1,
    "South GBR",
    ltmp_south,
    south_res;
    showlegend=true,
    legend_row=2,
    legend_col=2
)
linkyaxes!(ax1, ax2, ax3)

resize_to_layout!(f)

save_dir = OUT_DIR

mkpath(save_dir)

save("$(save_dir)/locs_reg.png", f)

f = taxa_cover_proportions(rs_raw.raw)

save("$(save_dir)/locs_taxa_cov.png", f)

f = taxa_population_proportions(rs_raw.raw)

save("$(save_dir)/locs_taxa_pop.png", f)

f = temporal_size_class_proportions(rs_raw.raw)

save("$(save_dir)/locs_size.png", f)

mkpath("$(save_dir)/loc_plots")

for i in 1:296
    location_comparison(rs_raw.raw, i, "$(save_dir)/loc_plots")
end
