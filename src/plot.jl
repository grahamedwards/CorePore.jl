Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" @eval using .CairoMakie
Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a"  @eval using .WGLMakie
Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" @eval using .WGLMakie

import CleanHistograms
export histograms, defaultlabels

include(download("https://raw.githubusercontent.com/grahamedwards/CleanHistograms.jl/main/src/CleanHistograms.jl"))

function table(x::NamedTuple)


end 

defaultlabels() = (onset = "Onset date (ka)", dfrz = "Freezing rate (m/yr)", dmlt = "Melting rate (m/yr)", sea2frz = "Basal freezing (‰, δ¹⁸O)", frz2mlt = "Basal melting (‰, δ¹⁸O)", flr = "Column depth (m)", basalCl = "Deep chloridity (g/kg)", basalO = "Deep δ¹⁸O (‰)")


function histograms(x::NamedTuple; f=Makie.Figure(size=(600,600)),  panels::NamedTuple = defaultlabels(), cols::Int=0, bins::Int=32, darkmode::Bool=false)
    n = length(panels)
    
    ks = keys(panels)

    pltclr = ifelse(darkmode,:white,:black)

    cols = ifelse(cols < 1, ceil(Int,sqrt(n)), cols)
    rows = ceil(Int,n/cols)

    @inbounds for i = 1:rows, j=1:cols
        ij = cols * (i-1) + j
        if ij < n
            k = ks[ij]
            @assert k ∈ keys(x) "Provided panel $k is not a key in chains provided in `x`"

            xk = x[k]

            ax =    Makie.Axis(f[i,j], xlabel=panels[k], 
                        bottomspinecolor=pltclr,xtickcolor=pltclr,xticklabelcolor=pltclr, xlabelcolor=pltclr,backgroundcolor=ifelse(darkmode,:transparent,:white),
                        xgridvisible=false,ygridvisible=false,yticklabelsvisible=false,yticksvisible=false,rightspinevisible=false,leftspinevisible=false,topspinevisible=false)
            h = CleanHistograms.cleanhist(xk,bins=bins)
            band!(ax,h.x,h.y,zero(h.y), color=(pltclr,0.1))
            lines!(ax,h.x,h.y, color=pltclr, linewidth=2,)
        end 
    end
    f
end



x=" "