using Test

include("../cover_construction.jl")

@testset "Exponential Average Area" begin
    lower_bound::Float64 = 0.0
    upper_bound::Float64 = 1.0

    lambda::Float64 = 1.0
    normalisation::Float64 = (1 - 1 / ℯ)

    @test average_area_expo(
        lambda, lower_bound, upper_bound
    ) == π / 4 * (2 - 5 / ℯ) / normalisation

    lower_bound = 1.0
    upper_bound = 3.0

    lambda = 2.0
    normalisation = (1 - exp(- lambda * upper_bound)) - (1 - exp(- lambda * lower_bound))

    @test average_area_expo(
        lambda, lower_bound, upper_bound
    ) == π / 2 * (5 / (4 * ℯ^2) - 25 / (4 * ℯ^6)) / normalisation
end

@testset "size class distribution" begin
    bounds::Vector{Float64} = [0.0, 1.0, 2.0, 3.0, 4.0]
    lambda = 1.0

    distribution = size_class_distribution(lambda, bounds)

    @test sum(distribution) ≈ 1.0
end

@testset "construct location cover" begin
    location_sample::Vector{Float64} = [
        0.85, # Location relative cover
        0.5,  # taxonomy weights
        0.5,
        0.5,
        0.5,
        0.5,
        1.0,  # Size distributions
        1.0,
        1.0,
        1.0,
        1.0,
    ]

    cover::Matrix{Float64} = zeros(Float64, 7, 5)
    bin_edges::Matrix{Float64} = [
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
    ]

    construct_location_cover!(
        cover,
        location_sample,
        bin_edges
    )
    @test sum(cover) ≈ 0.85
    @test all(cover[:, 1] .== cover[:, 2] .== cover[:, 3] .== cover[:, 4] .== cover[:, 5])
    @test all(sum(cover, dims=1) .≈ 0.85 / 5)
end
