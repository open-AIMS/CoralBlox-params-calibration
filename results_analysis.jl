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
    size_class_props = size_class_distribution(5.0, ADRIA.bin_edges()[1, :])
    for (idx, row_idx) in enumerate(ltmp_reefmod_idxs)
        loc_cov = rm_ltmp_taxa[2, :, idx] .* size_class_props' ./ sum(rm_ltmp_taxa[2, :, idx])
        tot_cov = raw_ltmp_reef_data[idx, findfirst(x->!ismissing(x), raw_ltmp_reef_data[idx, :])]
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

coral_param_fn = "C:/Users/dtan/data/coral_p_calib_fixed.dat"
coral_params = deserialize(coral_param_fn)

# Load values into scenario dataframe
scens = ADRIA.param_table(dom)
coral_param_values = coral_params[1:length(coral_p_names)]
scens[1, coral_p_names] = coral_param_values

# Extract and format location scale factors
scale_factors::Array{Float64, 3} = reshape(coral_params[length(coral_p_names)+1:end], (5, 4, 3))

rs_raw = ADRIA.run_model(dom, scens[1, :], scale_factors, dom_idxs)
