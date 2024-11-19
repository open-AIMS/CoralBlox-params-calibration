"""
Setup the parameter bounds and ordering for calibration.

Define standard gl0obal variables for which this ordering and naming can be accessed.
"""

"""
    target_param_names()
"""
function target_param_names()
    return [
        "linear_extension", "mb_rate", "fecundity", "dist_mean"
    ]
end

"""
    sc_fg_param_idxs(param_name::String, coral_df::DataFrame):mean_colony_diameter_m:Vector{Int64}

Return the indices of the parameter such that it returns them in the ordered used in ADRIA.

Makes it easier to insert bounds into vector of sample bounds from a ordered source file.

# Example

```julia
julia> ordered_idxs = sc_fg_param_idxs(param_name, coral_df)

julia> coral_df[ordered_idxs, :].fieldname
35-element Vector{Symbol}:
 :tabular_Acropora_1_1_linear_extension
 :tabular_Acropora_1_2_linear_extension
 :tabular_Acropora_1_3_linear_extension
...
 :corymbose_Acropora_2_1_linear_extension
...
 :corymbose_non_Acropora_3_1_linear_extension
...
 :small_massives_4_1_linear_extension
...
 :large_massives_5_1_linear_extension
...
```

"""
function sc_fg_param_idxs(param_name::String, coral_df::DataFrame)::Vector{Int64}
    idxs = zeros(Int64, 35)

    n_taxa, n_sizes = size(ADRIA.bin_widths())
    counter::Int64 = 1
    for grp in 1:n_taxa, sc in 1:n_sizes
        idxs[counter] = extract_param_group_idx(coral_df, "$(grp)_$(sc)_$(param_name)")[1]
        counter += 1
    end
    return idxs
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

"""
    set_bounds!(sample_bounds, coral_param_idxs, lower_bounds::Vector{Float64}, upper_bounds::Vector{Float64})::Nothing
"""
function set_bounds!(
    sample_bounds,
    coral_param_idxs,
    lower_bounds::Vector{Float64},
    upper_bounds::Vector{Float64}
)::Nothing
    if any(lower_bounds .> upper_bounds)
        msg = "At least one lower bound given is less then the corresponding upper bound"
        throw(ArgumentError(msg))
    end
    sample_bounds[coral_param_idxs] .= collect(zip(lower_bounds, upper_bounds))
    return nothing
end

function flatten_group_size(param_matrix::Matrix{Float64})::Vector{Float64}
    n_taxa, n_sizes = size(param_matrix)
    return reshape(permutedims(param_matrix, (2, 1)), (n_taxa * n_sizes,))
end

ecorrap_params = open_dataset(ECORRAP_PARAM_PATH)

# Construct linear extension parameter bounds
linear_extension_mean = ecorrap_params.filled_lin_ext_mean
linear_extension_mean[:, 7] .= 0.0

lin_ext_lb = linear_extension_mean.data[:, :] .* 0.9
lin_ext_ub = linear_extension_mean.data[:, :] .* 1.1

lin_ext_lb = flatten_group_size(lin_ext_lb)
lin_ext_ub = flatten_group_size(lin_ext_ub)

# Construct background mortality parameter bounds
mb_rate_mean = ecorrap_params.filled_mb_rate_mean

mb_rate_lb = mb_rate_mean.data[:, :] .* 0.9
mb_rate_ub = mb_rate_mean.data[:, :] .* 1.1

mb_rate_lb = flatten_group_size(mb_rate_lb)
mb_rate_ub = flatten_group_size(mb_rate_ub)

# Define parameter space to scan over
coral_params = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Coral)

# Extract just the target coral parameters
lin_ext_idx, mbrate_idx, fecundity_idx, dhw_tol_mean_idx =
    extract_param_group_idx.([coral_params], target_param_names())

coral_param_idx = vcat(
    lin_ext_idx, mbrate_idx, fecundity_idx, dhw_tol_mean_idx
)
coral_params = coral_params[sort(coral_param_idx), :]
coral_param_names = coral_params.fieldname

# Get updated parameter positions
lin_ext_idx, mbrate_idx, fecundity_idx, dhw_tol_mean_idx =
    extract_param_group_idx.([coral_params], target_param_names())

lin_ext_idx = sc_fg_param_idxs("linear_extension", coral_params)
mbrate_idx = sc_fg_param_idxs("mb_rate", coral_params)

sample_bounds = collect(zip(
    coral_params.lower_bound,
    coral_params.upper_bound
))

# Adjust bounds for fecundity and initial mean DHW tolerance
adjust_bounds!(sample_bounds, fecundity_idx, 0.5, 3.0)
adjust_bounds!(sample_bounds, dhw_tol_mean_idx, 0.8, 3.0)

# Set bounds for linear extension and mb rate from ecorrap data
set_bounds!(sample_bounds, lin_ext_idx, lin_ext_lb, lin_ext_ub)
set_bounds!(sample_bounds, mbrate_idx, mb_rate_lb, mb_rate_ub)

# Add parameters for location-specific scaling
n_groups = 5
n_size_classes = 7
n_limited_locs = length(limited_locations)
n_factors = 3  # growth, mortality, fecundity

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
        sampled_params[param_idxs[7]:param_idxs[8]], target_dom_idxs=loc_idxs
    )

    return new_dom, scen, growth_acc_params, scale_factors
end
