module PorewaterDiffusion
import DelimitedFiles 
using Random

export seawater, AND1B, AND2A, PorewaterProperty, SedimentColumn, LR04, constants, Proposal, update
include("parameters.jl")

export density, diffusionadvection, diffuseadvectcolumn!
include("diffuse.jl")

export porewaterhistory, porewaterhistory!
include("history.jl")

#export 
include("statistics.jl")

end
