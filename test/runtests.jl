using PorewaterDiffusion
using Test
using StableRNGs

macro silence(block)
    quote
        ose,oso = stderr,stdout
        redirect_stderr(devnull); redirect_stdout(devnull)
        x= $(esc(block))
        redirect_stderr(ose); redirect_stdout(oso)
        x
    end
end

mean(x) = sum(x)/length(x)

@testset "parameters" begin include("param.jl") end 

@testset "diffusion" begin include("diff.jl") end

@testset "histories" begin include("hist.jl") end

@testset "statistics" begin include("stats.jl") end

@testset "metropolis" begin include("metro.jl") end
