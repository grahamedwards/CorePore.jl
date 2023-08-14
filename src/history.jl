"""

```julia
porewaterhistory!(sc::SedimentColumn, p::Proposal, k::NamedTuple, climhist::NamedTuple, sw::NamedTuple, ka_dt::Int)
```

In-place version of [`porewaterhistory`](@ref), which takes every input as an arg (rather than some defaults as kwargs). It also requires you to provide `ka_dt` -- the number of diffusion timesteps in each thousand-year climate timestep.

see also: [`porewaterhistory`](@ref)

"""
function porewaterhistory!(sc::SedimentColumn, p::Proposal, k::NamedTuple, climhist::NamedTuple, sw::NamedTuple, ka_dt::Int)
    
    isd = searchsortedfirst(climhist.t, p.onset, rev=true)
    #isd = ifelse(isd<climhist.n, isd, climhist.n)
    @inbounds for t = isd:climhist.n
        benthic = climhist.x[t]

        @inbounds for j = 1:ka_dt

            Clo, Oo, ρ = boundaryconditions(sc.Cl.o[1], sc.O.o[1], benthic, p.sea2frz, p.frz2mlt, p.dmlt, p.dfrz, sw.Cl, sw.O, k.dz, k.dt)


            sc.Cl.o[1], sc.O.o[1], sc.rho.o[1] = Clo, Oo, ρ
            
            diffuseadvectcolumn!(sc,k)
        end
    end
    (; Cl = sc.Cl.p, d18O = sc.O.p)
end





"""

```julia
porewaterhistory(proposals [; k=constants(), climatehistory=LR04(), seawater=AND2A()])
```

Calculate the porewater advection-diffusion history of chlorinity and O-isotope-traced water in a sediment column described by properties in `k` (generated with [`constants`](@ref)) over a given [`climatehistory`](@ref) ([`LR04`](@ref) by default) and coretop `seawater` compositions.

The provided `proposals` (as custom type [`Proposal`](@ref)) describe the sensitivity and response of the system to climate fluctuations as recorded in `climatehistory`.

See [`diffuseadvectcolumn!`](@ref) for the underlying diffusion-advection transport calculations.

see also: [`porewaterhistory!`](@ref), [`Proposal`](@ref), [`constants`](@ref), [`LR04`](@ref), [`seawater`](@ref)

"""
function porewaterhistory(p::Proposal; k::NamedTuple=constants(), climatehistory::NamedTuple=LR04(), seawater::NamedTuple=AND2A())

    sc = SedimentColumn(k.nz,seawater...)
    x= porewaterhistory!(sc, p, k, climatehistory, seawater, dt_climatetimestep(climatehistory.t,k.dt))

    x
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