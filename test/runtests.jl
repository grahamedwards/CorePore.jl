using PorewaterDiffusion
using Test

@testset "custom types" begin 

@test Seawater(1,2).Cl == Seawater(1.,2.).Cl
@test AND1B().O == -0.33 && AND1B().Cl == 19.2657
@test AND2A().O == -1.0 && AND2A().Cl == 19.81655

@test length(PorewaterProperty(4).i) == 4
@test PorewaterProperty(2,1).o == PorewaterProperty([1.,1.],[1.,1.]).o

@test length(SedimentColumn(4).Cl.i) == 4
@test SedimentColumn(4, 1., 1., 1.).rho.o == SedimentColumn(PorewaterProperty(4,1.),PorewaterProperty(4,1.),PorewaterProperty(4,1.)).rho.o
end

@testset "diffusion" begin include("diffuse.jl") end
