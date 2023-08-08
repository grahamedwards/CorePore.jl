

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

`struct` to contain sediment column poperties at each node for a prior timestep `o` and current timestep `i`.

Constructor function returns an instance of `PorewaterProperty` with vectors of length `n`.  Optionally provide a value `x <: Number` to fill vectors with (otherwise values are undefined). 

"""
struct PorewaterProperty
    o::Vector{Float64}
    i::Vector{Float64}
end

PorewaterProperty(n::Int) = PorewaterProperty(Vector{Float64}(undef,n), Vector{Float64}(undef,n))

function PorewaterProperty(n::Int, x::Number)
    v=fill(float(x),n)
    PorewaterProperty(v,copy(v))
end





"""

    SedimentColumn(n::Int, [, Cl, O, rho])

`struct` to contain `PorewaterProperty`s for the porewater properties of [Cl⁻] (`Cl`), δ¹⁸O (`O`), and density (`rho`).

Constructor function returns an instance of `SedimentColumn` with `PorewaterProperty` vectors of length `n`.  Optionally provide values for `Cl`, `O`, and `rho` (otherwise values are undefined). 

see also: [`PorewaterProperty`](@ref)

"""
struct SedimentColumn
    Cl::PorewaterProperty
    O::PorewaterProperty
    rho::PorewaterProperty
end

SedimentColumn(n::Int) = SedimentColumn(PorewaterProperty(n), PorewaterProperty(n), PorewaterProperty(n))

SedimentColumn(n::Int, c::Number, o::Number, r::Number) = SedimentColumn(PorewaterProperty(n,c), PorewaterProperty(n,o), PorewaterProperty(n,r))