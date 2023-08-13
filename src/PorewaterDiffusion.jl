module PorewaterDiffusion

import DelimitedFiles 

export seawater, AND1B, AND2A, PorewaterProperty, SedimentColumn, LR04, constants, Proposal
include("parameters.jl")

export density, diffusionadvection, diffuseadvectcolumn!
include("diffuse.jl")

export porewaterhistory, porewaterhistory!
include("history.jl")
end
