module CorePore
import DelimitedFiles, Requires, Statistics
using Random

export Water, mcmurdoshelf, mcmurdosound, deepbonney, CoreData, andrill2a, andrill1b, PorewaterProperty, SedimentColumn, LR04, Constants, Proposal, update, ProposalPriors
include("parameters.jl")

export density, diffusionadvection, diffuseadvectcolumn!
include("diffuse.jl")

export porewaterhistory, porewaterhistory!
include("history.jl")

export loglikelihood, normpdf, means, medians
include("statistics.jl")

export porewatermetropolis
include("metropolis.jl")

function __init__()
    Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("plot.jl")
    Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("plot.jl") 
    Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" include("plot.jl") 
end
end
