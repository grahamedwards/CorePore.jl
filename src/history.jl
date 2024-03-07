"""

```julia
porewaterhistory!(sc::SedimentColumn, p::Proposal, k::Constants, climhist::NamedTuple, seawater::Water, ka_dt::Int)
```

In-place version of [`porewaterhistory`](@ref), which takes every input as an arg (rather than some defaults as kwargs). It also requires you to provide `ka_dt` -- the number of diffusion timesteps in each thousand-year climate timestep.

see also: [`porewaterhistory`](@ref)

"""
function porewaterhistory!(sc::SedimentColumn, p::Proposal, k::Constants, climhist::ClimateHistory, sw::Water, ka_dt::Int)
    #sc.Cl.p .= sc.Cl.o .= sw.Cl
    #sc.O.p .= sc.O.o .= sw.O
    equilibratecolumn!(sc,sw,water(p.basalCl,p.basalO),k.z,k.depth)

    isd = searchsortedfirst(climhist.t, p.onset, rev=true)
    #isd = ifelse(isd<climhist.n, isd, climhist.n)
    @inbounds for t = isd:climhist.n
        climate = climhist.x[t]

        @inbounds for j = 1:ka_dt

            Clo, Oo, ρ = boundaryconditions(sc.Cl.o[1], sc.O.o[1], climate, p.sea2frz, p.frz2mlt, p.dmlt, p.dfrz, sw.Cl, sw.O, k.dz, k.dt)

            sc.Cl.o[1], sc.O.o[1], sc.rho.o[1] = Clo, Oo, ρ
            
            diffuseadvectcolumn!(sc,k)
        end
    end
end





"""

```julia
porewaterhistory(p [; k=Constants(), climatehistory=LR04(), seawater=mcmurdosound()])
```

Calculate the porewater advection-diffusion history of chlorinity and O-isotope-traced water in a sediment column described by properties in `k` (::[`Constants`](@ref)) over a given [`ClimateHistory`](@ref) ([`LR04`](@ref) by default) and coretop `seawater` compositions.

[`Proposal`](@ref)) instance `p` describes the sensitivity and response of the system to climate fluctuations as recorded in `climatehistory`.

See [`diffuseadvectcolumn!`](@ref) for the underlying diffusion-advection transport calculations.

see also: [`porewaterhistory!`](@ref), [`Proposal`](@ref), [`Constants`](@ref), [`LR04`](@ref), [`water`](@ref)

"""
function porewaterhistory(p::Proposal; k::Constants=Constants(), climatehistory::ClimateHistory=LR04(), seawater::Water=mcmurdosound())

    sc = SedimentColumn(k.nz,seawater...)
    porewaterhistory!(sc, p, k, climatehistory, seawater, dt_climatetimestep(climatehistory.t,k.dt))

    (; Cl = sc.Cl.p, d18O = sc.O.p)
end





"""

    dt_climatetimestep(katime,dt)

A helper function to calculate the number of diffusion model timesteps `dt` (in years) in each timestep of the climate timescale `katime` (in kiloannum). Returns an integer. 

"""
function dt_climatetimestep(katime::AbstractRange,dt::Number) 
    y = abs(1000step(katime))
    @assert iszero(y%dt) "The timestep of katime (in years) is not divisible by dt"
    round(Int,y/dt)
end