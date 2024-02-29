module PorewaterDiffusion
import DelimitedFiles 
using Random

export Water, water, mcmurdoshelf, mcmurdosound, deepbonney, CoreData, andrill2a, PorewaterProperty, SedimentColumn, LR04, Constants, Proposal, update, getproposal
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
