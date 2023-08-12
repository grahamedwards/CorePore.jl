"""

```julia
porewaterhistory!(sc::SedimentColumn, p::Proposal, k::NamedTuple, climhist::NamedTuple, sw::Seawater, ka_dt::Int)
```

In-place version of [`porewaterhistory`](@ref), which takes every input as an arg (rather than some defaults as kwargs). It also requires you to provide `ka_dt` -- the number of diffusion timesteps in each thousand-year climate timestep.

see also: [`porewaterhistory`](@ref)

"""
function porewaterhistory!(sc::SedimentColumn, p::Proposal, k::NamedTuple, climhist::NamedTuple, sw::Seawater, ka_dt::Int)

    isd = searchsortedfirst(climhist.t, p.onset, rev=true)

    isd = ifelse(isd<climhist.n, isd, climhist.n)

    @inbounds for t = isd:climhist.n

        PorewaterDiffusion.boundaryconditions(sc.Cl.o[1], sc.O.o[1], climhist.x[t], p.sea2frz, p.frz2mlt, p.dmlt, p.dfrz, sw.Cl, sw.O, k.dz, k.dt)

        @inbounds for j = Base.OneTo(ka_dt)
            diffuseadvectcolumn!(sc,k)
        end
    end
end





"""

```julia
porewaterhistory(proposals [; k=constants(), climatehistory=LR04(), seawater=AND2A()])
```

Calculate the porewater advection-diffusion history of chlorinity and O-isotope-traced water in a sediment column described by properties in `k` (generated with [`constants`](@ref)) over a given [`climatehistory`](@ref) ([`LR04`](@ref) by default) and coretop `seawater` compositions.

The provided `proposals` (as custom type [`Proposal`](@ref)) describe the sensitivity and response of the system to climate fluctuations as recorded in `climatehistory`.

See [`diffuseadvectcolumn!`](@ref) for the underlying diffusion-advection transport calculations.

see also: [`porewaterhistory!`](@ref), [`Proposal`](@ref), [`constants`](@ref), [`LR04`](@ref), [`Seawater`](@ref)

"""
function porewaterhistory(p::Proposal; k::NamedTuple=constants(), climhist::NamedTuple=LR04(), sw::Seawater=AND2A())

    ka_dt =  round(Int,1000abs(step(climhist.t))/k.dt)
    sc = SedimentColumn(k.nz,sw.Cl, sw.O)
    porewaterhistory!(sc, p, k, climhist, sw, ka_dt)

    sc
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