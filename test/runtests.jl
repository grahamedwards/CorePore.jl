using PorewaterDiffusion
using Test

@testset "custom types" begin include("param.jl") end 

@testset "diffusion" begin include("diff.jl") end

@testset "histories" begin include("hist.jl") end
