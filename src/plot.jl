Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" @eval import .CairoMakie as Makie
Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a"  @eval import .GLMakie as Makie
Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" @eval import .WGLMakie as Makie

import Markdown
export table, histograms, defaultlabels, traces, sieveresults, profiledensity, profiletrace

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



"""

    sieveresults(chains; start=1, sieve=1, stop=auto, k=Constants(), climate=LR04(), seawater=mcmurdoshelf())

Sieve the chains of a [`porewatermetropolis`](@ref) run from `start` to `stop` (full chain length by default), including every `sieve`-th realization. Provide relevant Constants `k`, ClimateHistory `climate`, and `seawater` composition. 

Returns a NamedTuple of depth nodes `z`, and corresponding matrices of `Cl` and `O` where rows correspond to nodes in `z` and columns correspond to sieved results.

"""
function sieveresults(chains::NamedTuple; k::Constants=Constants(), climate::ClimateHistory=LR04(), seawater::Water = mcmurdoshelf(), sieve::Int=1, start::Int=1, stop::Int=0)

    chainlength = length(chains[1])
    @assert 0 < start < chainlength
    
    samples = start : sieve : ifelse(iszero(stop), chainlength, stop)
    @assert 0 <= last(samples) <= chainlength
    
    O = Matrix{Float64}(undef,k.nz, length(samples))
    Cl = similar(O)
    
    sc = SedimentColumn(k.nz,seawater...)
    ka_dt = CorePore.dt_climatetimestep(climate.t,k.dt)
    
    @inbounds for i = samples
        p = Proposal(chains.onset[i], chains.dfrz[i], chains.dmlt[i], chains.sea2frz[i], chains.frz2mlt[i], chains.flr[i], chains.basalCl[i], chains.basalO[i])
    
        porewaterhistory!(sc,p, k, climate, seawater, ka_dt)
    
        O[:,i] .= sc.O.p
        Cl[:,i] .= sc.Cl.p
    end
    
    return (; z= k.z, O, Cl)
end

"""

    coredensity(x, xz, prior, priorz; f=Figure(), bins=100, xlabel, ylabel, color)

Plot a heatmap of densities for sieved results `x` with corresponding depth nodes `xz`, and `prior` data at depths in `priorz`. Heatmap `bins`=100 by default, and you may optionally provide a custom `xlabel`, `ylabel`, and `color`map. See https://docs.makie.org/v0.21/explanations/colors for options. 

To add the plot as a panel in a prexisting figure, provide position information to `f` (e.g. `f=fig[1,2]`).

"""
function profiledenisty(x::Matrix{F}, xz::AbstractRange, prior::CorePore.MuSig, priorz::Vector; f=Makie.Figure(), bins::Int=100, xlabel::String="[Cl⁻] (g/kg)", ylabel::String="Depth below sea floor (m)", color=:Blues) where F<:AbstractFloat

    n = size(x,2)
    
    mn,mx = extrema(x)
    scooch = 0.1abs(mx-mn)
    br = range(mn-scooch,mx+scooch, bins+1)
    
    hists = Matrix{Float64}(undef, bins,length(xz))
    h = Vector{Int}(undef,bins)
    
    @inbounds for i = eachindex(xz)
        view(hists,:,i) .= CleanHistograms.quickhist!(h,view(x,i,:), br) / n
    end
    
    ax = Makie.Axis(f[1,1], yreversed=true,xlabel=xlabel, ylabel=ylabel)
    
    Makie.heatmap!(ax,br, xz, hists, lowclip=:transparent, colorrange=(0,1),colormap=color)
    
    Makie.errorbars!(ax,prior.mu, priorz, prior.sig, direction=:x, color=(:black), linewidth=1)

    Makie.scatter!(ax,prior.mu, priorz, color=:black, markersize=5)
    
    f
    end



"""

    coretrace(x, xz, prior, priorz; f=Figure(), bins=100, xlabel, ylabel, color)

Plot traces for sieved results `x` with corresponding depth nodes `xz` and`prior` data at depths in `priorz`. You may optionally provide a custom `xlabel`, `ylabel`, and trace `color`. See https://juliagraphics.github.io/Colors.jl/stable/namedcolors/ for options.

To add the plot as a panel in a prexisting figure, provide position information to `f` (e.g. `f=fig[1,2]`).

"""
function coretrace(x::Matrix{F}, xz::AbstractRange, prior::CorePore.MuSig, priorz::Vector; f=Makie.Figure(), xlabel::String="[Cl⁻] (g/kg)", ylabel::String="Depth below sea floor (m)", color=:tomato) where F<:AbstractFloat

    ax = Makie.Axis(f[1,1], yreversed=true,xlabel=xlabel, ylabel=ylabel)
    wt = 1/ size(x,2)

    @inbounds for i = axes(x,2)
        lines!(ax, view(x,:,i), xz, color=(color,wt))
    end

    Makie.errorbars!(ax,prior.mu, priorz, prior.sig, direction=:x, color=(:black), linewidth=1)
    Makie.scatter!(ax,prior.mu, priorz, color=:black, markersize=5)

    f
end