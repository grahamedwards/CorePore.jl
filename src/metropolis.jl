"""

    stopwatch(i, n, t)

Convenience function for [`porewatermetropolis`] that returns a String reporting the progress at step `i` for total steps `n` with start time `t` (in s since the epoch).

"""
function stopwatch(i::Integer,n::Integer,t::Number)
    pd = 10 * i ÷ n
    bar = "■"^pd * "□"^(10-pd)
    tt = round((time() - t) / 60,digits=2)
    string("0% |", bar,"| 100%  ||  step: $i / $n  ||  time: $tt m")
end

"""

    PorewaterDiffusion.strictpriors(p::Proposal, record_max_age::Number, climatelimits::Tuple{Number,Number}, depth::Number)

Evalute strict constraints on priors that will automatically reject a proposal with...
- Onset date beyond the climate record timespan (in ka) (`record_max_age`)
- Nonphysical subglacial thresholds -- melting at lower benthic δ¹⁸O than freezing or values exceeding the record extrema (`climatelimits`).
- Freezing rate is non-zero and less than the upperbound of observed freezing rates (0.002 m/yr)
- Annual melting rate is non-zero and must be no more than that observed at the Thwaites grounding line (<10 m/yr, [Davis+ 2023](https://www.nature.com/articles/s41586-022-05586-0)).
- Diffusive porewater column (`p.flr`) is between 0 and the model sediment column `depth`.
- Basal [Cl⁻] compositions < 0 g/kg or in excess of 200 g/kg (the water+halite peritectic is ~163 g/kg Cl, assuming charge balance with NaCl)
-  δ¹⁸O compositions less than dome-like values (~ -56 ‰)

"""
function strictpriors(p::Proposal, record_max_age::Number, climatelimits::Tuple{Number,Number}, depth::Number)
    
    x=true
    
    x &= 0 < p.onset <= record_max_age
    x &= climatelimits[1] < p.sea2frz < p.frz2mlt < climatelimits[2]
    x &= 0 < p.dfrz <= 0.002 # ≤ 0.4dt/dz to prevent an error in a log-calculation in `PorewaterDiffusion.boundaryconditions`
    x &= 0 < p.dmlt <= 10.
    x &= 0 < p.flr <= depth
    x &= 0 < p.basalCl <= 200
    x &= -56 < p.basalO

    x
end



"""

    PorewaterDiffusion.proposaljump(p::Proposal, σ::Proposal; f=proposals, rng::AbstractRNG)

Add a random jump to a randomly selected field of `p` with a corresponding normal jumping distribution defined by the corresponding field in `σ`. The possible fields may be specified by providing a Tuple of Symbols `f`, and a specific RNG seed may be provided.

Note that `:dmlt`, `:dfrz`, and `:basalCl` are drawn from a lognormal jumping distribution (where `σ` is in log-space.)

"""
function proposaljump(p::Proposal, j::Proposal; rng::AbstractRNG=Xoshiro(), f::Tuple{Vararg{Symbol}}=proposals)

    jumpname = rand(rng,f)
    logdist = jumpname ∈ (:dmlt, :dfrz, :basalCl)
    jump = j[jumpname] * randn(rng)
    x = p[jumpname]
    x = logdist ? log(x) : x
    x += jump
    x = ifelse(logdist, exp(x),x)
    (update(p, jumpname, x) , jumpname , abs(jump) )
end



"""

```julia
porewatermetropolis...
```
Not tested, yet...
"""
function porewatermetropolis(p::Proposal, jumpsigma::Proposal, prior::CoreData; burnin::Int=0, chainsteps::Int=100, k::Constants=Constants(), seawater::Water=mcmurdosound(), explore::Tuple{Vararg{Symbol}}=proposals, climate::ClimateHistory=LR04(), rng::AbstractRNG=Random.Xoshiro())

    scalejump=2.4

    record_max_age = first(climate.t)
    climate_limits = extrema(climate.x)

    ϕ = p # make a new proposal from the original.

    chains = Matrix{Float64}(undef, length(proposals), chainsteps)
    lldist = Vector{Float64}(undef, chainsteps)
    acceptance = falses(chainsteps)

    sc = SedimentColumn(k.nz,seawater...)
    ka_dt = PorewaterDiffusion.dt_climatetimestep(climate.t,k.dt)
    
    porewaterhistory!(sc, ϕ, k, climate, seawater, ka_dt)
    
    llCl, llO = loglikelihood(prior.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p), loglikelihood(prior.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)
    ll = llCl + llO
    
    clock = time()
    burnupdate, chainupdate = div.((ifelse(iszero(burnin),1,burnin), ifelse(iszero(chainsteps),1,chainsteps)),20,RoundUp)
    println("Beginning sequence...\n  $burnin burn-in iterations \n  $chainsteps recorded iterations\n ------------ \n\n " )
    flush(stdout)

    burninacceptance=0
    
    @inbounds for i=Base.OneTo(burnin)


        ϕ, jumpname, jump = proposaljump(p, jumpsigma, f=explore, rng=rng)
        if strictpriors(ϕ, record_max_age, climate_limits, k.depth)

            porewaterhistory!(sc, ϕ, k, climate, seawater, ka_dt)

            llCl, llO = loglikelihood(prior.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p), loglikelihood(prior.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)
            llϕ = llCl + llO
        else
            llϕ=-Inf
        end

        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            jumpsigma = update(jumpsigma, jumpname, jump * scalejump) # update jumpsigma
            p, ll = ϕ, llϕ  # update proposal and log-likelihood
            burninacceptance=+1        
        end

        if iszero(i % burnupdate) # Update progress
            println("Burn-In --- ", stopwatch(i,burnin,clock))
            flush(stdout)
        end
    end

    
    println("\n\n$burnin burn-in steps complete. ℓ = $ll, acceptance rate= $(100burninacceptance÷ifelse(iszero(burnin),1,burnin)) %.\n\nCurrent guess: $p\nJumps = $jumpsigma\n"); flush(stdout)


    @inbounds for i = 1:chainsteps

        ϕ, jumpname, jump = proposaljump(p, jumpsigma, f=explore, rng=rng)
        if strictpriors(ϕ, record_max_age, climate_limits, k.depth)

            porewaterhistory!(sc, ϕ, k, climate, seawater, ka_dt)

            llCl, llO = loglikelihood(prior.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p), loglikelihood(prior.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)
            llϕ = llCl + llO
        else
            llϕ=-Inf
        end
        
        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            jumpsigma = update(jumpsigma, jumpname, jump * scalejump) # update jumpsigma
            p, ll = ϕ, llϕ  # update proposal and log-likelihood
            acceptance[i] = true   
        end
        

        chains[:,i] .= (p...,)
        lldist[i] = ll

        if iszero(i % chainupdate) # Update progress 
            println("Main Chain --- ", stopwatch(i,chainsteps,clock))
            flush(stdout)
        end

    end
    outnames = (explore...,:ll, :accept)
    outvalues = ((chains[i,:] for i in axes(chains,1))..., lldist, acceptance)
    NamedTuple{outnames}(outvalues)
end


