Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" @eval import .CairoMakie as Makie
Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a"  @eval import .WGLMakie as Makie
Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" @eval import .WGLMakie as Makie

import Markdown
export table, histograms, defaultlabels, traces

include(download("https://raw.githubusercontent.com/grahamedwards/CleanHistograms.jl/main/src/CleanHistograms.jl"))


"""

    table(x)

Print a styled Markdown table in the REPL of the summary statistics in `x`, generated with [`means`](@ref) or [`medians`](@ref).

"""
function table(x::NamedTuple)
    k = keys(x)
    h = keys(x[1])

    t = if length(h)==1
        "||m|\n|:--|:--|"
    elseif length(h)==2
        "|| μ | σ |\n|:--|:--|:--|"
    elseif length(h)==3
        "|| M | – | + |\n|:--|:--|:--|:--|"
    end

    @inbounds for i = keys(x)
        l = "|$i"
        xi = x[i]
        @inbounds for j = eachindex(xi)
            l *= "|$(xi[j])"
        end
        t *= "\n$l|"
    end
    Markdown.parse(t)
end 


"""

    defaultlabels()

Returns default labels corresponding to fields in [`Proposal`](@ref).

"""
defaultlabels() = (onset = "Onset date (ka)", dfrz = "Freezing rate (m/yr)", dmlt = "Melting rate (m/yr)", sea2frz = "Basal freezing (‰, δ¹⁸O)", frz2mlt = "Basal melting (‰, δ¹⁸O)", flr = "Column depth (m)", basalCl = "Deep chloridity (g/kg)", basalO = "Deep δ¹⁸O (‰)")


"""

    histograms(x; f=Makie.Figure(), panels=defaultlabels(), bins=32, cols=auto, darkmode=false)

Returns panels of histograms for NamedTuple of posterior chains `x` over `col` columns. Designate specific fields of `x` to plot in panels as a `NamedTuple` with names in `x` paired to labels as `String`s. Plots all panels by default.
Optionally invert colorscheme with darkmode=`true`.

"""
function histograms(x::NamedTuple; f=Makie.Figure(size=(600,600)),  panels::NamedTuple = defaultlabels(), cols::Int=0, bins::Int=32, darkmode::Bool=false)

    n, ks = length(panels), keys(panels)

    pltclr = ifelse(darkmode,:white,:black)

    cols = ifelse(cols < 1, ceil(Int,sqrt(n)), cols)
    rows = ceil(Int,n/cols)

    @inbounds for i = 1:rows, j=1:cols
        ij = cols * (i-1) + j
        if ij < n
            k = ks[ij]
            @assert k ∈ keys(x) "Provided panel $k is not a key in chains provided in `x`"

            xk = x[k]

            ax = Makie.Axis(f[i,j], xlabel=panels[k], 
                    bottomspinecolor=pltclr,xtickcolor=pltclr,xticklabelcolor=pltclr, xlabelcolor=pltclr,backgroundcolor=ifelse(darkmode,:transparent,:white),
                    xgridvisible=false,ygridvisible=false,yticklabelsvisible=false,yticksvisible=false,rightspinevisible=false,leftspinevisible=false,topspinevisible=false)
            h = CleanHistograms.cleanhist(xk,bins=bins)
            Makie.band!(ax,h.x,h.y,zero(h.y), color=(pltclr,0.1))
            Makie.lines!(ax,h.x,h.y, color=pltclr, linewidth=2,)
        end 
    end
    f
end



"""

   traces(x; f=Makie.Figure(), domain, panels=defaultlabels(), cols=auto, darkmode=false)

Returns panels of trace plots for NamedTuple of posterior chains `x` over `col` columns. Designate specific fields of `x` to plot in panels as a `NamedTuple` with names in `x` paired to labels as `String`s. Plots all panels by default. Optionally provide a UnitRange of Markov chain steps to plot with `domain`. Optionally invert colorscheme with darkmode=`true`.

"""
function traces(x; f=Makie.Figure(size=(800,400)), domain::UnitRange= 0:0, panels=defaultlabels(), cols::Int=0, darkmode=false)
    n, ks = length(panels), keys(panels)

    d = ifelse(length(domain)>1, domain, 1:length(x[1]))

    pltclr = ifelse(darkmode,:white,:black)

    cols = ifelse(cols < 1, ceil(Int,sqrt(n)), cols)
    rows = ceil(Int,n/cols)

    @inbounds for i = 1:rows, j=1:cols
        ij = cols * (i-1) + j
        if ij < n
            k = ks[ij]
            @assert k ∈ keys(x) "Provided panel $k is not a key in chains provided in `x`"

            xk = x[k]

            ax = Makie.Axis(f[i,j], ylabel=panels[k], 
                    leftspinecolor=pltclr,ytickcolor=pltclr,yticklabelcolor=pltclr, ylabelcolor=pltclr,backgroundcolor=ifelse(darkmode,:transparent,:white),
                    xgridvisible=false,ygridvisible=false,xticklabelsvisible=false,xticksvisible=false,rightspinevisible=false,bottomspinevisible=false,topspinevisible=false)
            Makie.lines!(ax,xk[d],color=pltclr)
        end 
    end
    f
end