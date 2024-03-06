using PorewaterDiffusion
using Test
using StableRNGs



@testset "parameters" begin include("param.jl") end 

@testset "diffusion" begin include("diff.jl") end

@testset "histories" begin include("hist.jl") end

@testset "statistics" begin include("stats.jl") end

@testset "metropolis" begin include("metro.jl") end
