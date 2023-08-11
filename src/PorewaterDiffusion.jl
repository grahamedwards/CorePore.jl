module PorewaterDiffusion

import DelimitedFiles 

export Seawater, AND1B, AND2A, PorewaterProperty, SedimentColumn, LR04, constants
include("parameters.jl")


include("diffuse.jl")

end
