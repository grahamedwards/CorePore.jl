

"""

    water(Cl, O)

Custom `struct` containing `Float64` values of chlorinity `Cl` and δ¹⁸O `O`. Built-in splatting functionality.

see also: [`mcmurdoshelf`](@ref), [`mcmurdosound`](@ref)

"""
@kwdef struct Water
    Cl::Float64
    O::Float64
end

# Add splatting functionality
function Base.iterate(@nospecialize(x::Water),i::Int=1)
    y = ifelse(i==1, x.Cl, nothing)
    y = ifelse(i==2, x.O, y)
    ifelse(isnothing(y), nothing, (y,i+1))
end

"""
    mcmurdoshelf()

Generate a [`Water`](@ref) instance with coretop porewater compositions from the McMurdo ice shelf. Pairs with core ANDRILL-1B.

see also: [`andrill1b`](@ref)

"""
mcmurdoshelf() = Water(19.2657,-0.33)

"""
    mcmurdosound()

Generate a [`Water`](@ref) instance with southern McMurdo Sound water compositions. Pair with core ANDRILL-2A.

see also:  [`andrill2a`](@ref)

"""
mcmurdosound() = Water(19.81655,-1.0)

"""
    deepbonney()

    Generate a [`Water`](@ref) instance with >30 m Lake Bonney water compositions. Use for sediment column basal boundary condition. Chloridity (143.333 g/L  ÷ 1.2 kg/L ≈ 119 g/kg) from Angino+ 1963 ([doi:10.1086/626879](https://doi.org/10.1086/626879)) and δ¹⁸O from Matsubaya+ 1979 ([doi:10.1016/0016-7037(79)90042-5](https://doi.org/10.1016/0016-7037(79)90042-5)).
    

"""
deepbonney() = Water(119.44416666666667,-25.2)




@kwdef struct MuSig
    mu::Vector{Float64}
    sig::Vector{Float64}
end

"""

    CoreData(z, mCl, σCl, mO, σO)

Returns an  `CoreData` struct with sediment core data formatted for [`porewatermetropolis`](@ref) — `(; z, Cl = (mu, sig), O = (mu, sig))`. Inputs are Vectors for values of sample depths `z` (meters below sea floor), measured (mean) chlorinity `mCl` and 1σ uncertainty `σCl`, and measured (mean) δ¹⁸O `mO` and 1σ uncertainties `σO`. 

Vectors must all be of the same length. If Cl and δ¹⁸O sampling is not 1:1, use `NaN` or anything that is not a subtype of Number (e.g. `missing`, `nothing`). While `NaN` is used internally, `CoreData` does this conversion for you.

---

    CoreData(z, m, σ, measurement)

Same as above, but for a core with only chlorinity or δ¹⁸O data (the other Vectors are empty to return null values in [`loglikelihood`](@ref) calculations). Provide sample depths `z`, measured means `m`, 1σ uncertainties `σ`, and the `measurement` as a symbol, e.g. `:Cl` or `:O`. 

---

`CorePore.jl` comes loaded with convenience functions to generate data for [`andrill2a`](@ref) from Tracy+ 2010 ([doi:10.1130/G30849.1](https://doi.org/10.1130/G30849.1)) and [`andrill1b`](@ref) from Pompilio+ 2007 (https://digitalcommons.unl.edu/andrillrespub/37). 

"""
struct CoreData
    z::Vector{Float64}
    Cl::MuSig
    O::MuSig
end


function CoreData(z::Vector{<:Number}, mCl::AbstractVector, σCl::AbstractVector, μO::AbstractVector, σO::AbstractVector)
    @assert length(z) == length(mCl) == length(σCl) == length(μO) == length(σO)

    @inbounds for i= eachindex(z)
        mCl[i] = ifelse(typeof(mCl[i])<:Number, mCl[i], NaN)
        σCl[i] = ifelse(typeof(σCl[i])<:Number, σCl[i], NaN)
        μO[i] = ifelse(typeof(μO[i])<:Number, μO[i], NaN)
        σO[i] = ifelse(typeof(σO[i])<:Number, σO[i], NaN)
    end

    CoreData(float.(z), MuSig(float.(mCl), float.(σCl)), MuSig(float.(μO), float.(σO)))
end
function CoreData(z::Vector{<:Number}, m::AbstractVector, σ::AbstractVector, s::Symbol)
    
    @assert length(z) == length(m) == length(σ)
    
    chlorine, oxygen = (:Cl, :CL, :chlorine, :chlorinity), (:O, :oxygen, :water, :d18O, :δ18O, :δ¹⁸O, :d¹⁸O) 
    
    @assert s ∈ (chlorine..., oxygen...)

    @inbounds for i= eachindex(z)
        m[i] = ifelse(typeof(m[i])<:Number, m[i], NaN)
        σ[i] = ifelse(typeof(σ[i])<:Number, σ[i], NaN)
    end
    x = MuSig(float.(m), float.(σ))
    xx = Vector{eltype(x.mu)}(undef,0)
    xx = MuSig(xx, xx)

    s ∈ chlorine ? CoreData(float.(z), x, xx) : CoreData(float.(z), xx, x)
end


"""

    andrill2a

Generate a [`CoreData`](@ref) instance with data from core ANDRILL-2A (from Frank+ 2010, [doi:10.1130/G30849.1](https://doi.org/10.1130/G30849.1)). 

see also: [`CoreData`](@ref)

"""
function andrill2a()
    z = [9.67, 30.09, 37.41, 43.72, 51.3, 57.21, 62.66, 73.15, 81.03, 92.97, 116.22, 155.76, 235.66, 293.3, 336.18, 353.53, 545.01, 619.35, 779.69, 809.84, 963.44]

    Clm = (35.45/1000) .* [654, 576, 612, 659, 693, 692, 691, 821, 722, 740, 804, 1117, 1974, 2100, 2303, 2253, 2771, 2728, 2895, 2722, 3091]

    Cls = .03Clm # 3-5 % reported reproducibility. We choose the lower bound given the ±2% precision reported in Pompilio+ 2007.

    Om = [-1.3, -2.7, -5.6, -8.1, -9.8, -10.0, -10.6, -10.3, missing, -10.9, -5.2, -8.5, -9.7, -10.6, -10.2, missing, -9.3, missing, missing, missing, missing]

    Os = fill(0.1, length(Om)) # based on UCSC isotope lab.

    CoreData(z, Clm,Cls, Om, Os)
end



"""

    andrill1b

Generate a [`CoreData`](@ref) instance with data from core ANDRILL-1B (from Pompilio+ 2007 (https://digitalcommons.unl.edu/andrillrespub/37)). 

see also: [`CoreData`](@ref)

"""
function andrill1b()
    z = [9.95, 20.55, 30.82, 39.36, 54.06, 69.36, 79.59, 96.75, 112.51, 138.76, 189.775, 223.955, 251.45, 284.3, 375.19, 410.25, 477.8, 519.96, 572.685, 610.105, 649.575, 668.205, 798.055, 848.985, 901.775, 936.78, 967.58, 1076.72]

    Clm = (35.45/1000) .* [475.1914646, 645.7566202, 680.8842872, 712.4045889, 792.2377042, 849.2664369, 810.2374237, 811.5119651, 859.4921896, 797.2260501, 868.1681727, 823.9844697, 816.0152075, 828.4942182, 935.8452392, 901.5187553, 918.5586697, 889.6388733, 838.8985985, 759.6680245, 743.8599926, 658.1991176, 624.5028023, 593.2801276, 563.1605987, 552.8358149, 568.6385208, 485.6941693]

    Cls = .02Clm # ±2% precision reported

    Om = [missing, missing, missing, -9.29, -8.94, -7.56, -7.94, -8.17, -3.03, -1.99, -8.24, -8.48, -7.66, -7.37, -7.35, -7.09, -6.29, -6.21, -5.92, -3.65, missing, missing, missing, missing, -1.74, -1.62, missing, missing]
    
    Os = fill(0.1, length(Om)) # based on UCSC isotope lab.

    CoreData(z, Clm,Cls, Om, Os)
end


"""

    PorewaterProperty(n::Int, [, x])

`struct` to contain sediment column poperties at each node for the previous timestep `o` and present timestep `p`.

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

`struct` to contain `PorewaterProperty`s for the porewater properties of chloridity (`Cl`), δ¹⁸O (`O`), and density (`rho`).

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

Generate a [`ClimateHistory`](@ref) instance containing the Liesiecki & Raymo 2004 benthic stack ([data](https://lorraine-lisiecki.com/stack.html), [publication](https://doi.org/10.1029/2004PA001071)) interpolated for 1 ka timesteps and going forward in time from 5.32 Ma. 

| field | description |
| :---- | :---------- | 
| `t`   | time (ka)  |
| `x`   | benthic δ¹⁸O (‰) |
| `n`   | timesteps in `t` |

"""
function LR04()
    a = DelimitedFiles.readdlm(string(@__DIR__,"/../data/LR04-interpolated-1ka.csv"),',')
    @assert a[2,1] - a[1,1] ≈ 1
    t = a[end,1] : -1 : a[1,1]
    x = reverse(a[:,2])
    ClimateHistory(t, x, length(t))
end




"""

    ClimateHistory

Custom struct to contain a climate record `x` over timescale `t`, where `t` contains `n` timesteps.
    
The type generated by [`LR04`](@ref), which is used for (perhaps unnecessary) type stability and code readability.

"""
@kwdef struct ClimateHistory
    t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int}
    x::Vector{Float64}
    n::Int
end





"""

```julia
Constants( ; k, dt, dz, depth )
```

Returns a `Constants` struct containing constants and coefficients used in diffusion history calculations. The inputs and their default values are listed in the table below. From these it calculates a few convenient variables---the product `dt * dz` (`dtdz`), a Range of node depths `z` (in m), the number of nodes `nz`, and the `penultimate_node`---as well as temperature-dependent diffusion coefficients used in `diffusionadvection`, using the temperature-depth paramterization of [Morin+ 2010](https://doi.org/10.1130/GES00512.1). These are returned as Vectors of length `nz`, where each cell corresponds to a depth node in `z`. The two coefficients are `k1` and `k2`, both of which are follwed with `cl` or `w` to denote Cl⁻ or water, respectively: `k1cl`, `k2cl`, `k1w`, and `k2w`.

| field | description | default |
| :---- | :---------- | ----: |
| `k`   | hydraulic conductivity of sediment column | 0.1 m / yr
| `dt` | timestep | 10 yrs |
| `dz` | node spacing | 5 m |
| `depth` | sediment column depth | 2000 m |

"""
struct Constants

    k::Float64
    dz::Float64
    dt::Float64
    depth::Float64
    dtdz::Float64
    z::AbstractRange{Float64}
    nz::Int
    penultimate_node::Int
    k1cl::Vector{Float64}
    k2cl::Vector{Float64}
    k1w::Vector{Float64}
    k2w::Vector{Float64}
end

function Constants(; k::Number=0.1, dt::Number=10., dz::Number=5.,  depth::Number=2000.)

    depth = ifelse( iszero(depth % dz), depth, depth - depth % dz)

    k, dz, dt = float.((k, dz, dt))

    z = 0.0 : dz : depth
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

    Constants(k, dz, dt, depth, dt*dz, z, nz, nz-1, k1cl, k2cl, k1w, k2w)
end 


"""

    Proposal

Struct containing porewater parameters (described below). All inputs must be of type Number (converts to Float64). Built-in functionality includes splatting (`...`) and [`CorePore.fastsplat`](@ref).

| field | description | units | default |
| :---- | :---------- | :---- | :------ |
|`onset`| onset of model | ka | 5320 |
| `dfrz`| freezing rate | m/yr | 0.0001 | 
| `dmlt`| melting rate | m/yr | 0.001 |
| `sea2frz` | Benthic δ¹⁸O threshold for subglacial freezing | ‰ | 3.5 |
| `frz2mlt` | Benthic δ¹⁸O threshold for subglacial melting | ‰ | 4.2 |
| `flr` | depth of diffusion-dominated porewater column | m | 1000. |
| `basalCl` | chloridity at base of diffusion-dominated column | g/kg | 119.4 |
| `basalO` | δ¹⁸O at base of diffusion-dominated column | ‰ | -25.2 |

see also: [`update`](@ref)

"""
@kwdef struct Proposal
    onset::Float64 = 5320.
    dfrz::Float64 = 1e-4
    dmlt::Float64 = 1e-3
    sea2frz::Float64 = 3.5
    frz2mlt::Float64 = 4.2
    flr::Float64 = 1000.
    basalCl::Float64 = 119.4
    basalO::Float64 = -25.2
end

# Add splatting functionality for Proposal
function Base.iterate(@nospecialize(x::Proposal),i::Int=1)
    y = ifelse(i==1, x.onset, nothing)
    y = ifelse(i==2, x.dfrz, y)
    y = ifelse(i==3, x.dmlt, y)
    y = ifelse(i==4, x.sea2frz, y)
    y = ifelse(i==5, x.frz2mlt, y)
    y = ifelse(i==6, x.flr, y)
    y = ifelse(i==7, x.basalCl, y)
    y = ifelse(i==8, x.basalO, y)
    ifelse(isnothing(y), nothing, (y,i+1))
end


"""
    update(p::Proposal, f::Symbol, x::Number)

Update field `f` of [`Proposal`](@ref) instance `p` with value `x`.

"""
function update(x::Proposal, f::Symbol,v::Number)
    @assert f ∈ fieldnames(Proposal)
    v = float(v)
    
    Proposal(
            ifelse(f==:onset, v, x.onset),
            ifelse(f==:dfrz, v, x.dfrz),
            ifelse(f==:dmlt, v, x.dmlt),
            ifelse(f==:sea2frz, v, x.sea2frz),
            ifelse(f==:frz2mlt, v, x.frz2mlt),
            ifelse(f==:flr, v, x.flr),
            ifelse(f==:basalCl, v, x.basalCl),
            ifelse(f==:basalO, v, x.basalO)
    )
end

"""
    fastsplat(x::Proposal)

Returns a tuple of the contents of `x`. Avoids type inherent instability of `Base.iterate` for fast splatting.
"""
fastsplat(x::Proposal) = (x.onset, x.dfrz, x.dmlt, x.sea2frz, x.frz2mlt, x.flr, x.basalCl, x.basalO)



"""

    ProposalPriors(climate, k; onset, dfrz, fmlt, climatelimits, flr, basalCl, basalO)


Custom `struct` to hold the bounds of parameters in Proposal as `Tuple`s of (`minimum value`, `maximum value`). The constructor function takes a [`ClimateHistory`](@ref) instance `climate` and [`Constants`](@ref) instance `k`, as well as any customizations to default values.

---

## Fields

See [`Proposal`](@ref) for descriptions and units of shared parameters.

| field | default | explanation |
| :---- | :------ | :---------- |
|`onset`|  (`0`, `first(climate.t)`) | Onset date must fall within climate record timespan.
| `dfrz`| (`0`, `0.002`) | Freezing rate >0 and less than the upperbound of observed freezing rates (0.002 m/yr).¹ |
| `dmlt` | (`0`, `10`) | Annual melting rate >0 and less than that observed at the Thwaites grounding line.² |
| `climatelimits` | `extrema(climate.x)` | `sea2frz` and `frz2mlt` must lie within the observed climate record values. |
| `flr` | (`0`, `k.depth`) | Diffusive porewater column (`p.flr`) has non-zero depth within model domain. |
| `basalCl` | (`0`, `200`) | Basal [Cl⁻] compositions within 0-200 g/kg. (³) |
| `basalO` | (`-56`, `∞`) | δ¹⁸O exceed dome-like values. |
|||

¹ Must be  ≤ 0.4dt/dz to prevent an error in a log-calculation in `CorePore.boundaryconditions`

² <10 m/yr -- Davis+ 2023, https://www.nature.com/articles/s41586-022-05586-0

³ The water+halite peritectic is ~163 g/kg Cl, assuming charge balance with NaCl.

"""
struct ProposalPriors
    onset::Tuple{Float64,Float64}
    dfrz::Tuple{Float64,Float64}
    dmlt::Tuple{Float64,Float64}
    climatelimits::Tuple{Float64,Float64} 
    flr::Tuple{Float64,Float64}
    basalCl::Tuple{Float64,Float64}
    basalO::Tuple{Float64,Float64} 
end

function ProposalPriors(climate::ClimateHistory, k::Constants; 
    onset::Tuple{Number,Number}=(NaN,NaN),
    dfrz::Tuple{Number,Number} = (0.,0.002),
    dmlt::Tuple{Number,Number} = (0.,10.),
     climatelimits::Tuple{Number,Number} = (NaN,NaN),
     flr::Tuple{Number,Number}  = (NaN,NaN), 
     basalCl::Tuple{Number,Number} = (0.,200.), 
     basalO::Tuple{Number,Number} = (-56.,Inf))

    onset = ifelse(isnan(sum(onset)), (0,first(climate.t)), onset)
    climatelimits = ifelse(isnan(sum(climatelimits)), extrema(climate.x), climatelimits)
    flr = ifelse(isnan(sum(flr)), (0,k.depth), flr)
  
    ProposalPriors(float.(onset), float.(dfrz), float.(dmlt), float.(climatelimits), float.(flr), float.(basalCl), float.(basalO))
end