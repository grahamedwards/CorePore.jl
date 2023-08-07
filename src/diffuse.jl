
"""

```julia
function boundaryconditions(Cl, d18O, x, ocean2freeze,freeze2melt, meltrate, freezerate, seawater, dz, dt)
```
Calculates sediment surface boundary condition for δ¹⁸O (`d18O`) and [Cl⁻] (`Cl`), based on the thermodynamic state described by the current current benthic δ¹⁸O value `x` and the threshold values corresponding to subglacial freezing `ocean2freeze` and subglacial melting `freeze2melt`.

For melting or freezing states, calculates boundary condition from the assumed `meltrate`, `freezerate`, timestep `dt`, length-step `dt`, and composition of `seawater` (NamedTuple).

"""
function boundaryconditions(Cl, d18O, x, ocean2freeze,freeze2melt, meltrate, freezerate, seawater, dz, dt)
    
    @assert ocean2freeze <= freeze2melt
    
    if x < ocean2freeze # low δ18O -> warm -> seawater
        Cl = seawater.cCl
        d18O = seawater.d18O

    elseif x > freeze2melt # high δ18O -> cold -> warm-based
        ϕdz = dz*0.4
        melt = meltrate * dt
        
        Cl *= ϕdz / (ϕdz + melt)
        d18O *= ϕdz * melt / (ϕdz+melt)
    
    else # mid δ18O -> mid temps -> cold-based
        ϕdz = 0.4dz # scaled for porosity
        frz = freezerate*dt
    
        Cl *= ϕdz / (ϕdz - frz) 
        d18O +=  1.59 * log(1 - (frz / ϕdz)) # simplified from eqn 2 of Toyota et al. (2017)
    end
    (Cl, d18O)
end
