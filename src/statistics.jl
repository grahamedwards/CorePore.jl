"""
```julia 
normpdf(x, μ, σ)
```
Calculate the probability density of a normal distribution with mean `μ` and standard deviation `σ` at the value `x`.
"""
normpdf(μ::Number,σ::Number,x::Number) = exp(-(x-μ)*(x-μ) / (2*σ*σ)) / (σ*sqrt(2*π))





"""
```julia 
normll(x, μ, σ)
```
Calculate the (relative) log-likelihood that an observation `x` was drawn from the normal distribution with mean `μ` and standard deviation `σ`. 

Note: This function excludes a constant that will not vary for different vlaues of `x` to speed up calculation in `metropolis`. To calculate the absolute log-likelihood, take the `log` of `normpdf`
"""
normll(μ::Number,σ::Number,x::Number) = -(x-μ)*(x-μ) / (2*σ*σ)





"""

    linterp(x, y₂, y₁, x₁, Δx)

Estimate the value of y corresponding to `x` given known coordinate (`x₁`, `y₁`), value `y₂`, and `Δx`= x₂ - x₁. Note that if `x`==`x₁`, y= `y₁`. 

"""
linterp(x,x1,dx,y2,y1) = y1 + (x-x1) * (y2-y1)/dx






"""
    loglikelihood( zₒ, μ, σ, zₘ, m ) 

Calculate the (relative) log-likelihood that model values in `m` simulated at depth nodes in `zₘ` were drawn from normally distributed observations with sample depths, mean values, and 1σ uncertainties in corresponding indices of the Vectors `z`, `μₒ`, and `σ`. In instances where sample depths do not correspond to depth nodes, the simulated value is linearly interpolated between bounding depth nodes in `zₘ`.

see also: [`normll`](@ref), [`linterp`](@ref)

"""
function loglikelihood(zo::T, muo::T, sigo::T, zm::AbstractRange, m::T) where T<:Vector{Float64}
    ll=zero(eltype(m))
    dzm = step(zm)

    @inbounds @simd for i = eachindex(zo)
        mi = searchsortedfirst(zm, zo[i])
        xmi = linterp(zo[i], zm[mi], dzm, m[mi+1], m[mi])
        ll += normll(muo[i], sigo[i], xmi)
    end
    ll
end





"""

    strictpriors(p, record_max_age, benthiclimits)

Evalute strict constraints on priors that will automatically reject a proposal:
- Onset date beyond the climate record timespan (in ka)
- Negative freezing or melting rates.
- Nonphysical subglacial thresholds -- melting at lower benthic δ¹⁸O than freezing or values exceeding the record extrema.

"""
function strictpriors(p::Proposal, record_max_age::Number, benthiclimits::Tuple)
    x=true
    x *= p.onset < record_max_age
    x *= 0 < p.onset

    x *= p.dfrz > 0 
    x *= p.dmlt > 0 
    
    x *= p.sea2frz < p.frz2mlt
    x *= benthiclimits[1] < p.sea2frz
    x *= p.frz2mlt < benthiclimits[2]

    x
end



function proposaljump(p::Proposal, j::Proposal; f::NTuple=fieldnames(Proposal), rng::AbstractRNG=Random.Xoshiro())

    jumpname = rand(rng,f)
    jump = getproperty(j,jumpname) * randn(rng)
    jumped = getproperty(p,jumpname) + jump

    (update(p,jumpname,jumped) , jumpname , jump )
end





"""

```julia
porewatermetropolis...
```
Not tested, not exported... yet...
"""
function porewatermetropolis(p::Proposal, jumpsize::Proposal, prior::NamedTuple; burnin::Int=0, chainsteps::Int=100, explore=fieldnames(Proposal), k::NamedTuple=constants(), seawater::NamedTuple=AND2A(), benthic::NamedTuple=LR04(), scalefactor=2.9, rng::AbstractRNG=Random.Xoshiro())

    record_max_age = first(benthic.t)
    benthic_limits = extrema(benthic.x)

    ϕ = p # make a new proposal from the original.

    chains = Matrix{eltype(p)}(undef, chainsteps, length(fieldnames(Proposal)))

    sc = SedimentColumn(k.nz,seawater...)
    ka_dt = PorewaterDiffusion.dt_climatetimestep(ch.t,k.dt)
    
    porewaterhistory!(sc, ϕ, k, ch,sw, ka_dt)

    ll= loglikelihood(prior.Cl.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p) + loglikelihood(prior.O.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)

    @inbounds for i=Base.OneTo(burnin)

        ϕ, jumpname,jump = proposaljump(p, jumpsize)
        if strictpriors(ϕ, record_max_age, benthic_limits)

            porewaterhistory!(sc, ϕ, k, ch, sw, ka_dt)

            llϕ = loglikelihood(prior.Cl.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p) + loglikelihood(prior.O.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)
        else
            llϕ=-Inf
        end

        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            jumpsize = update(jumpsize,jumpname,abs(jump)*stepfactor) # update jumpsize
            p = ϕ  # update proposal
            ll = llϕ # Record new log likelihood              
        end
    end

    @inbounds for i=Base.OneTo(chainsteps)

        ϕ, jumpname,jump = proposaljump(p, jumpsize)
        if strictpriors(ϕ, record_max_age, benthic_limits)

            porewaterhistory!(sc, ϕ, k, ch, sw, ka_dt)

            llϕ = loglikelihood(prior.Cl.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p) + loglikelihood(prior.O.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)
        else
            llϕ=-Inf
        end

        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            jumpsize = update(jumpsize,jumpname,abs(jump)*stepfactor) # update jumpsize
            p = ϕ  # update proposal
            ll = llϕ # Record new log likelihood              
        end

        chains[i,:] = p.onset, p.dfrz, p.dmlt, p.sez2frz, p.frz2mlt
    end
    NamedTuple{fieldnames(Proposal)}(chains[:,i] for i in axes(chains,2))
end