include("common.jl")
include("cover_construction.jl")
include("1_setup.jl")

# ----- CONSTRUCT PARAMETER NAMES IN CORRECT ORDER AS USED -----

sample_bounds = []

taxa_names = ["tabular_Acropora", "corymbose_Acropora", "corymbose_non_Acropora", "small_massives", "large_massives"]
szs = 1:7

coral_p_names = []

# Create growth bounds
size_widths = ADRIA.bin_widths()
for (t_idx, taxa) in enumerate(taxa_names)
    for s in szs
        push!(coral_p_names, Symbol(taxa * "_" * string(t_idx) * "_" * string(s) * "_" * "linear_extension"))
    end
end

# add mortality bounds
for (t_idx, taxa) in enumerate(taxa_names)
    for s in szs
        push!(coral_p_names, Symbol(taxa * "_" * string(t_idx) * "_" * string(s) * "_" * "mb_rate"))
    end
end


# add fecundity bounds
for (t_idx, taxa) in enumerate(taxa_names)
    for s in 4:7
        push!(coral_p_names, Symbol(taxa * "_" * string(t_idx) * "_" * string(s) * "_" * "fecundity"))
    end
end

# ----- LOAD Initial Coral Cover -----
function insert_init_loc_cover!(
    dom;
    raw_ltmp_reef_data=raw_ltmp_reef_data,
    rm_ltmp_taxa=rm_ltmp_taxa,
    ltmp_reefmod_idxs=ltmp_reefmod_idxs
)::Nothing
    size_class_props = size_class_distribution(0.5, ADRIA.bin_edges()[1, :])
    for (idx, row_idx) in enumerate(ltmp_reefmod_idxs)
        loc_cov = rm_ltmp_taxa[2, :, idx] .* size_class_props' ./ sum(rm_ltmp_taxa[2, :, idx])
        tot_cov = raw_ltmp_reef_data[idx, findfirst(x->!ismissing(x), raw_ltmp_reef_data[idx, :])] ./ dom.loc_data.k[row_idx]
        dom.init_coral_cover[:, row_idx] .= reshape(permutedims(loc_cov, (2, 1)), (35,)) .* tot_cov
    end
    return nothing
end

init_cover_fn = "C:/Users/dtan/data/init_cover.dat"
init_cover = deserialize(init_cover_fn)

construct_cover!(
    dom, init_cover, location_classification.consecutive_classification
)

insert_init_loc_cover!(dom)


# ----- LOAD CALIBRATED RESULTS -----

coral_param_fn = "C:/Users/dtan/repos/coral_p_calib_fixed.dat"
coral_params = deserialize(coral_param_fn)

# Load values into scenario dataframe
scens = ADRIA.param_table(dom)
coral_param_values = coral_params[1:length(coral_p_names)]
scens[1, coral_p_names] = coral_param_values

# Extract and format location scale factors
scale_factors::Array{Float64, 3} = reshape(coral_params[length(coral_p_names)+1:end], (5, 4, 3))

rs_raw = ADRIA.run_model(dom, scens[1, :], scale_factors, dom_idxs)

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

prefix = "corrected_growth_min_dhw"

save_dir = "Outputs/$(prefix)"

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
