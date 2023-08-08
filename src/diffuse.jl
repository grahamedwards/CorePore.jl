# Diffusion functions and support
# boundaryconditions
# 



"""

```julia
function boundaryconditions(Cl, d18O, x, ocean2freeze,freeze2melt, meltrate, freezerate, seawater, dz, dt)
```
Calculates sediment surface boundary condition for δ¹⁸O (`d18O`) and [Cl⁻] (`Cl`), based on the thermodynamic state described by the current current benthic δ¹⁸O value `x` and the threshold values corresponding to subglacial freezing `ocean2freeze` and subglacial melting `freeze2melt`.

For melting or freezing states, calculates boundary condition from the assumed `meltrate`, `freezerate`, timestep `dt`, length-step `dt`, and composition of seawater `Clsw` and `d18Osw`.

"""
function boundaryconditions(Cl, d18O, x, ocean2freeze,freeze2melt, meltrate, freezerate, Clsw, d18Osw, dz, dt)
    
    @assert ocean2freeze <= freeze2melt
    
    if x < ocean2freeze # low δ18O -> warm -> seawater
        Cl = Clsw
        d18O = d18Osw

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



"""

    diffusionadvection(x,above,below,k1,k2,v,dt,dz)

Calculate the property of a node in a vertical profile given the combined effects of diffusion and advection. Returns the property given initial values for the node `x`, the overlying node `above`, the underlying node `below`, diffusion coefficients `k1` and `k2`, vertical advection velocity `v`, timestep `dt`, and lengthstep `dz`. Alternatively provide the product of `v * dt * dz` for a minor speed-up.

"""
diffusionadvection(x,above,below,k1,k2,v,dt,dz) = x + k1 * (above - 2x + below) + k2 * (above - x) - (v*dt*dz) * (x - above)

diffusionadvection(x,above,below,k1,k2,vdtdz) = x + k1 * (above - 2x + below) + k2 * (above - x) - vdtdz * (x - above)


"""

    density(chlorinity)

Calculates the density of a water parcel with `chlorinity` in units g/m³ (rather than g/m³ for convenience with `velocity`)

see also: [`velocity`](@ref)

"""
density(chlorinity) = (chlorinity * 0.0018) + 1





"""

    velocity(x, above, k)

Calculate the velocity (m/yr) at a node with density `x`, given the density of the node `above`, and the hydraulic conductivity `k` (m/yr).

see also: [`density`](@ref)

"""
velocity(x, above, k) = k * (above - x) / x

