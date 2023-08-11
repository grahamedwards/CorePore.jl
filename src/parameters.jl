

"""

    Seawater(Cl, O)

Custom struct to hold chlorinity `Cl` and δ¹⁸O `O` values (necessarily `Float64`s). 

see also: [`AND1B`](@ref), [`AND2A`](@ref)

"""

struct Seawater
    Cl::Float64
    O::Float64
end

Seawater(cl::Number,o::Number) = Seawater(float(cl),float(o))

"""
    AND1B()

Generate a `Seawater` instance with coretop values of ANDRILL-1B.

see also: [`Seawater`](@ref)

"""
AND1B() = Seawater(19.2657,-0.33)

"""
    AND2A()

Generate a `Seawater` instance with coretop values of ANDRILL-2A.

see also: [`Seawater`](@ref)

"""
AND2A() = Seawater(19.81655,-1.0)





"""

    PorewaterProperty(n::Int, [, x])

`struct` to contain sediment column poperties at each node for a prior timestep `o` and present timestep `p`.

Constructor function returns an instance of `PorewaterProperty` with vectors of length `n`.  Optionally provide a value `x <: Number` to fill vectors with (otherwise values are undefined). 

"""
struct PorewaterProperty
    o::Vector{Float64}
    p::Vector{Float64}
end

PorewaterProperty(n::Int) = PorewaterProperty(Vector{Float64}(undef,n), Vector{Float64}(undef,n))

function PorewaterProperty(n::Int, x::Number)
    v=fill(float(x),n)
    PorewaterProperty(v,copy(v))
end





"""

    SedimentColumn(n::Int, [, Cl, O])

`struct` to contain `PorewaterProperty`s for the porewater properties of [Cl⁻] (`Cl`), δ¹⁸O (`O`), and density (`rho`).

Constructor function returns an instance of `SedimentColumn` with `PorewaterProperty` vectors of length `n`.  Optionally provide values for `Cl`, `O`, and `rho` (otherwise values are undefined). 

see also: [`PorewaterProperty`](@ref), [`density`](@ref)

"""
struct SedimentColumn
    Cl::PorewaterProperty
    O::PorewaterProperty
    rho::PorewaterProperty
end

SedimentColumn(n::Int) = SedimentColumn(PorewaterProperty(n), PorewaterProperty(n), PorewaterProperty(n))

SedimentColumn(n::Int, c::Number, o::Number) = SedimentColumn(PorewaterProperty(n,c), PorewaterProperty(n,o), PorewaterProperty(n,density(c)))





"""

    LR04()

Load a NamedTuple containing the Liesiecki & Raymo 2004 benthic stack ([data](https://lorraine-lisiecki.com/stack.html), [publication](https://doi.org/10.1029/2004PA001071)) interpolated for 1 ka timesteps and going forward in time from 5.32 Ma. 

| field | description |
| :---- | :---------- | 
| `t`   | time (ka)  |
| `x`   | benthic δ¹⁸O (‰)|

"""
function LR04()
    a = DelimitedFiles.readdlm(string(@__DIR__,"/../data/LR04-interpolated-1ka.csv"),',')
    @assert a[2,1] - a[1,1] ≈ 1
    t = a[end,1] : -1 : a[1,1]
    x = reverse(a[:,2])
    (; t, x)
end





"""

```julia
    constants( ; k, dt, dz, depth )
````

Returns a NamedTuple of constants and coefficients used in diffusion history calculations. The input constants and their default values are listed in the table below. From these it calculates a few convenience variables: the number of nodes `nz`, a range of `interiornodes`,  and the product `dtdz`. 

Using the constants, it calculates the number of nodes  as well as temperature-dependent diffusion coefficients used in `diffusionadvection`, using the temperature-depth paramterization of [Morin+ 2010](https://doi.org/10.1130/GES00512.1). These are returned as Vectors of length `nz`, where each cell corresponds to a depth node. The two coefficients are `k1` and `k2`, both of which are follwed with `cl` or `w` to denote Cl⁻ or water, respectively: `k1cl`, `k2cl`, `k1w`, and `k2w`.

| field | description | default |
| :---- | :---------- | ----: |
| `k`   | hydraulic conductivity of sediment column | 0.1 m / yr
| `dt` | timestep | 10 yrs |
| `dz` | node spacing | 5 m |
| `depth` | sediment column depth | 2000 m |

"""
function constants(; k::Number=0.1, dt::Number=10., dz::Number=5.,  depth::Number=2000.)

    depth = ifelse( iszero(depth % dz), depth, depth - depth % dz)

    k, dz, dt = float.((k, dz, dt))

    z = 0:dz:depth
    depth=last(z)
    nz = length(z)

    T = @. (0.0767 * z ) + 270.75  # Morin et al. 2010 || in K  (270.75  = 273.15 - 2.4)
    κCl = @. exp(3.8817 - (2.2854e3 / T)) # m² / yr
    κwater = @. exp(4.2049 - ( 2.2699e3 / T)) # m² / yr

    k1cl = κCl .* (dt / (dz * dz))
    k1w = κwater  .* (dt / (dz * dz))

    k2cl = Vector{eltype(T)}(undef,nz)
    k2w = similar(k2cl)
    k2w[1] = k2cl[1] = zero(eltype(T))

    @inbounds for i = 2:nz
        x = dt / dz
        k2cl[i] = (κCl[i] - κCl[i-1]) * x
        k2w[i] = (κwater[i] - κwater[i-1]) * x 
    end

    (; k, dz, dt, dtdz = dt*dz, depth, nz, k1cl, k2cl, k1w, k2w, interiornodes=2:nz-1)
end 