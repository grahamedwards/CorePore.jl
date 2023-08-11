module PorewaterDiffusion

import DelimitedFiles 

export Seawater, AND1B, AND2A, PorewaterProperty, SedimentColumn, LR04, constants, Proposal
include("parameters.jl")

export density, diffusionadvection, diffuseadvectcolumn!
include("diffuse.jl")

export porewater 
include("time.jl")
end
