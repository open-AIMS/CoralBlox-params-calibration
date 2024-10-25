using Test
using Distributions, Statistics

include("../src/common/cover_construction.jl")

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
        0.25,  # Size distributions
        0.25,
        0.25,
        0.25,
        0.25,
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

@testset "random cover construction" begin
    n_tests::Int64 = 100
    cover::Matrix{Float64} = zeros(Float64, 7, 5)
    bin_edges::Matrix{Float64} = [
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
        0.0 0.05 0.075 0.1 0.2  0.4 1.0 1.5;
    ]
    for _ in 1:n_tests
        location_sample = rand(Uniform(0,1), 11)

        construct_location_cover!(
            cover,
            location_sample,
            bin_edges
        )
        @test sum(cover) ≈ location_sample[1]
    end

end
