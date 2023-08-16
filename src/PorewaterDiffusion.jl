module PorewaterDiffusion
import DelimitedFiles 
using Random

export seawater, mcmurdoshelf, mcmurdosound, coredata, andrill2a, PorewaterProperty, SedimentColumn, LR04, constants, Proposal, update
include("parameters.jl")

export density, diffusionadvection, diffuseadvectcolumn!
include("diffuse.jl")

export porewaterhistory, porewaterhistory!
include("history.jl")

export loglikelihood, normpdf
include("statistics.jl")

export porewatermetropolis
include("metropolis.jl")

end
