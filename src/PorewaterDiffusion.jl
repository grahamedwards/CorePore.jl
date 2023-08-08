module PorewaterDiffusion

export Seawater, AND1B, AND2A, PorewaterProperty, SedimentColumn
include("parameters.jl")


include("diffuse.jl")

end
