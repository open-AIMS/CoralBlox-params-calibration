"""
Setup the parameter bounds and ordering for calibration.

Define standard gl0obal variables for which this ordering and naming can be accessed.
"""

"""
    target_param_names()
"""
function target_param_names()
    return [
        "linear_extension", "mb_rate", "mean_colony_diameter_m", "fecundity", "dist_mean"
    ]
end

"""
    adjust_bounds!(sample_bounds, coral_param_idx, scale_lb::Float64, scale_ub::Float64)::Nothing

Scale the given parameter bounds by the given scale factors inplace.
"""
function adjust_bounds!(
    sample_bounds, coral_param_idx, scale_lb::Float64, scale_ub::Float64
)::Nothing
    extended_lb = first.(sample_bounds[coral_param_idx]) .* scale_lb
    extended_ub = last.(sample_bounds[coral_param_idx]) .* scale_ub
    sample_bounds[coral_param_idx] .= collect(zip(extended_lb, extended_ub))
    return nothing
end

# Define parameter space to scan over
coral_params = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Coral)

# Extract just the target coral parameters
lin_ext_idx, mbrate_idx, coldiam_idx, fecundity_idx, dhw_tol_mean_idx =
    extract_param_group_idx.([coral_params], target_param_names())

coral_param_idx = vcat(
    lin_ext_idx, mbrate_idx, coldiam_idx, fecundity_idx, dhw_tol_mean_idx
)
coral_params = coral_params[sort(coral_param_idx), :]
coral_param_names = coral_params.fieldname

# Get updated parameter positions
lin_ext_idx, mbrate_idx, coldiam_idx, fecundity_idx, dhw_tol_mean_idx =
    extract_param_group_idx.([coral_params], target_param_names())

sample_bounds = collect(zip(
    coral_params.lower_bound,
    coral_params.upper_bound
))

# Adjust bounds for linear extensions, fecundity and initial mean DHW tolerance
adjust_bounds!(sample_bounds, lin_ext_idx, 0.25, 1.8)
adjust_bounds!(sample_bounds, fecundity_idx, 0.5, 3.0)
adjust_bounds!(sample_bounds, dhw_tol_mean_idx, 0.8, 3.0)

# Add parameters for location-specific scaling
n_groups = 5
n_size_classes = 7
n_limited_locs = length(limited_locations)
n_factors = 3  # growth, mortality, fecundity

# Adjust bounds for mortality rate
# Size specific changes for each group
for grp in 1:n_groups, sz in 1:n_size_classes
    group_idx = extract_param_group_idx(coral_params, "$(grp)_$(sz)_mb_rate")

    if sz < 4
        sample_bounds[group_idx] .= fill((0.01, 0.6), length(group_idx))
    else
        ub = grp > 3 ? 0.15 : 0.4
        sample_bounds[group_idx] .= fill((0.01, ub), length(group_idx))
    end
end

location_coef = fill((0.3, 1.5), n_groups * n_limited_locs * n_factors)

coral_start_idx = 1
coral_end_idx = length(sample_bounds)

loc_coef_start_idx = coral_end_idx + 1
append!(sample_bounds, location_coef)
loc_coef_end_idx = length(sample_bounds)

# Location-based growth scaling
growth_acc_start_idx = loc_coef_end_idx + 1
for _ in 1:length(TARGET_DOM_IDXS)
    push!(sample_bounds, (-30.0, -15.0))  # steepness
    push!(sample_bounds, (0.0, 2.0))  # height
    push!(sample_bounds, (0.0, 0.3))  # midpoint
end
growth_acc_end_idx = length(sample_bounds)

sc_dist_bounds = fill((0.25, 30.0), n_limited_locs)

sc_dist_start_idx = growth_acc_end_idx + 1
append!(sample_bounds, sc_dist_bounds)
sc_dist_end_idx = length(sample_bounds)

global PARAM_IDXS = [
    coral_start_idx, coral_end_idx,
    loc_coef_start_idx, loc_coef_end_idx,
    growth_acc_start_idx, growth_acc_end_idx,
    sc_dist_start_idx, sc_dist_end_idx
]

global CORAL_PARAM_NAMES = coral_params.fieldname

# Utility functions for domain and parameter calibration setup

"""
Reshape the growth acceleration parameters from a vector a matrix of shape [n_parameters ⋅ n_locs]
"""
function reshape_growth_accel_parameters(
    params::Vector{Float64};
    n_locs=length(TARGET_DOM_IDXS)
)::Matrix{Float64}
    return reshape(params, (3, n_locs))
end

"""
    insert_init_loc_cover!(dom, lambdas; raw_ltmp_reef_data=raw_ltmp_reef_data, rm_ltmp_taxa=rm_ltmp_taxa, target_dom_idxs=TARGET_DOM_IDXS )::Nothing

Recalculate initial cover of target locations to match photogrammetry coral composition and
manta tow total cover. Use the given lambda calculation
"""
function insert_init_loc_cover!(
    dom,
    lambdas;
    raw_ltmp_reef_data=raw_ltmp_reef_data,
    rm_ltmp_taxa=rm_ltmp_taxa,
    target_dom_idxs=TARGET_DOM_IDXS
)::Nothing
    for (idx, row_idx) in enumerate(TARGET_DOM_IDXS)
        size_class_props = size_class_distribution(lambdas[idx], ADRIA.bin_edges()[1, :])
        loc_cov =
            rm_ltmp_taxa[2, :, idx] .* size_class_props' ./ sum(rm_ltmp_taxa[2, :, idx])
        tot_cov =
            raw_ltmp_reef_data[
                idx, findfirst(x -> !ismissing(x), raw_ltmp_reef_data[idx, :])
            ] ./ dom.loc_data.k[row_idx]
        dom.init_coral_cover[:, row_idx] .=
            reshape(permutedims(loc_cov, (2, 1)), (35,)) .* tot_cov
    end

    return nothing
end

"""
    setup_run(dom::Domain, sampled_params::Vector{Float64}; param_names::Vector{Symbol}=CORAL_PARAM_NAMES, param_idxs=PARAM_IDXS, loc_idxs=target_dom_idxs )::Nothing

Insert coral parameters into dataframe, reconstruct size class distribution of target
locations and extract scale factors and growth acceleration parameters.
"""
function setup_run(
    dom::Domain,
    sampled_params::Vector{Float64};
    param_names::Vector{Symbol}=CORAL_PARAM_NAMES,
    param_idxs=PARAM_IDXS,
    loc_idxs=TARGET_DOM_IDXS
)::Tuple{Domain, DataFrame, Matrix{Float64}, Array{Float64, 3}}
    new_dom = deepcopy(dom)

    scen = ADRIA.param_table(new_dom)
    coral_param_values = sampled_params[param_idxs[1]:param_idxs[2]]
    scen[1, param_names] = coral_param_values

    growth_acc_params = reshape_growth_accel_parameters(
        sampled_params[param_idxs[5]:param_idxs[6]]
    )

    scale_factors::Array{Float64,3} = reshape(
        sampled_params[param_idxs[3]:param_idxs[4]], (5, 4, 3)
    )

    insert_init_loc_cover!(
        new_dom,
        sampled_params[param_idxs[7]:param_idxs[8]],
        target_dom_idxs=loc_idxs
    )

    return new_dom, scen, growth_acc_params, scale_factors
end
