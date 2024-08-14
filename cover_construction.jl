using ADRIA

"""
Average area of coral given a truncated exponential distribution. This is the expected value of
the π r^2
"""
function average_area_expo(lambda::Float64, lower_bound::Float64, upper_bound::Float64)::Float64
    lb_exp_term::Float64 = exp(-lambda * lower_bound)
    ub_exp_term::Float64 = exp(-lambda * upper_bound)

    lb::Float64  = lower_bound * lambda
    lb2::Float64 = (lb)^2

    ub::Float64  = upper_bound * lambda
    ub2::Float64 = (ub)^2

    normalisation::Float64 = (1 - ub_exp_term) - (1 - lb_exp_term)

    return π / 4 * (ub_exp_term / lambda^2 * ( -ub2 + 2 * ( -ub - 1)) -
                lb_exp_term / lambda^2 * ( -lb2 + 2 * ( -lb - 1))) / normalisation
end

function size_class_distribution(lambda::Float64, bin_edges::AbstractVector{Float64})::Vector{Float64}
    # Calculate average area assuming a truncated exponential distribution = E[πr^2]
    area_props::Vector{Float64} = average_area_expo.(lambda, bin_edges[1:end-1], bin_edges[2:end])
    area_props = area_props ./ sum(area_props)
    return area_props
end

"""
    construct_location_cover!(preallocated, location_sample, bin_edges)::Array{Float64, 2}

Location Sample element description

1        - Location relative habitable cover
[2 - 6]  - Taxonomy cover weights
[7 - 11] - Taxonomy size exponential parametrisation

# Arguments
- `preallocated` : An abstract Matrix of size [n_sizes ⋅ n_taxa]
- `location_sample` : An abstract vector containing the embedded from the sampler
- `bin_edges` : An abstract vector defining the diameter bounds of the taxa
"""
function construct_location_cover!(
    preallocated::AbstractMatrix{Float64},
    location_sample::AbstractVector{Float64},
    bin_edges::AbstractMatrix{Float64}
)::Nothing
    n_taxa::Int64 = 5
    # Calculate relative cover for each taxonomy and reuse location sample memory
    @views location_sample[2:(1 + n_taxa)] .*= (
        location_sample[1] ./ sum(location_sample[2:(1 + n_taxa)])
    )

    # Calculate size class weightings for each taxonomy
    for taxa in 1:n_taxa
        preallocated[:, taxa] .= size_class_distribution(location_sample[6 + taxa], bin_edges[taxa, :])
    end

    # Multiply size class weightings and taxonomy relative cover to create cover for
    # location
    preallocated .*= location_sample[2:(1 + n_taxa)]'
    return nothing
end

"""
First Float describes relative habitable cover. Next 5 Floats describe Taxonomy weightings,
next 5 descibe size class exponential paramterisation.
"""
function construct_cover!(dom::Domain, vec_sample::Vector{Float64}, location_types::Vector{Int64})::Nothing
    n_location_types = maximum(location_types)
    temporary_cover::Matrix{Float64} = zeros(Float64, 7, 5)
    bin_edges::Matrix{Float64} = ADRIA.bin_edges()

    location_mask::BitVector = BitVector([true for _ in 1:3806])

    stride::Int64 = 11
    for loc_type in 1:n_location_types
        @views construct_location_cover!(
            temporary_cover,
            vec_sample[(1 +kstride * (loc_type - 1)):(stride * loc_type)],
            bin_edges
        )
        location_mask .= location_types .== loc_type
        dom.init_coral_cover[:, location_mask] .= reshape(temporary_cover, (35,))
    end
    return nothing
end
