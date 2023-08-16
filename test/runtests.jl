using PorewaterDiffusion
using Test
using StableRNGs
using Suppressor: @suppress



@testset "custom types" begin include("param.jl") end 

@testset "diffusion" begin include("diff.jl") end

@testset "histories" begin include("hist.jl") end

@testset "statistics" begin include("stats.jl") end
