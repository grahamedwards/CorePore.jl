

"""

    strictpriors(p, record_max_age, climatelimits)

Evalute strict constraints on priors that will automatically reject a proposal:
- Onset date beyond the climate record timespan (in ka)
- Negative freezing or melting rates.
- Nonphysical subglacial thresholds -- melting at lower benthic δ¹⁸O than freezing or values exceeding the record extrema.

"""
function strictpriors(p::Proposal, record_max_age::Number, climatelimits::Tuple)
    x=true
    x *= p.onset <= record_max_age
    x *= 0 < p.onset

    x *= p.dfrz > 0 
    x *= p.dmlt > 0 
    
    x *= p.sea2frz < p.frz2mlt
    x *= climatelimits[1] < p.sea2frz
    x *= p.frz2mlt < climatelimits[2]

    x
end



function proposaljump(p::Proposal, j::Proposal; rng::AbstractRNG, f::Tuple=fieldnames(Proposal))

    jumpname::Symbol = rand(rng,f)
    jump = getproposal(j,jumpname) * randn(rng)
    jumped = getproposal(p,jumpname) + jump

    (update(p,jumpname,jumped) , jumpname , jump )
end





"""

    stopwatch(i, n, clock)

Convenience function for [`porewatermetropolis`] that returns a String reporting the progress at step `i` for total steps `n` with start time `t` (in s since the epoch).

"""
function stopwatch(i::Integer,n::Integer,clock::Number)
    pd = 10 * i ÷ n
    bar = "■"^pd * "□"^(10-pd)
    t = round((time() - clock) / 60,digits=2)
    string("0% |", bar,"| 100%  ||  step: $i / $n  ||  time: $t m")
end





"""

```julia
porewatermetropolis...
```
Not tested, yet...
"""
function porewatermetropolis(p::Proposal, jumpsize::Proposal, prior::CoreData; burnin::Int=0, chainsteps::Int=100, k::Constants=Constants(), seawater::Seawater=mcmurdosound(), explore::Tuple=fieldnames(Proposal),climate::ClimateHistory=LR04(), scalejump=1.8, rng::AbstractRNG=Random.Xoshiro())

    record_max_age = first(climate.t)
    climate_limits = extrema(climate.x)

    ϕ = p # make a new proposal from the original.

    chains = Matrix{Float64}(undef, length(fieldnames(Proposal)), chainsteps)
    lldist = Vector{Float64}(undef, chainsteps)
    acceptance = falses(chainsteps)

    sc = SedimentColumn(k.nz,seawater...)
    ka_dt = PorewaterDiffusion.dt_climatetimestep(climate.t,k.dt)
    
    porewaterhistory!(sc, ϕ, k, climate, seawater, ka_dt)

    ll= loglikelihood(prior.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p) + loglikelihood(prior.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)

    clock, burnupdate, chainupdate = time(), div(burnin,20), div(chainsteps,20)
    println("Beginning sequence...\n  $burnin burn-in iterations \n  $chainsteps recorded iterations\n ------------ \n\n " )
    flush(stdout)

    burninacceptance=0
    
    @inbounds for i=Base.OneTo(burnin)

        ϕ, jumpname, jump = proposaljump(p, jumpsize, rng=rng)
        if strictpriors(ϕ, record_max_age, climate_limits)

            porewaterhistory!(sc, ϕ, k, climate, seawater, ka_dt)

            llCl = loglikelihood(prior.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p) 
            llO = loglikelihood(prior.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)
            llϕ = llCl + llO
        else
            llϕ=-Inf
        end

        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            jumpsize = update(jumpsize,jumpname,abs(jump)*scalejump) # update jumpsize
            p = ϕ  # update proposal
            ll = llϕ # Record new log likelihood
            burninacceptance=+1              
        end

        if iszero(i % burnupdate) # Update progress
            println("Burn-In --- ", stopwatch(i,burnin,clock))
            flush(stdout)
        end
    end

    println("\n\n$burnin burn-in steps complete. ℓ = $ll, acceptance rate= $(100burninacceptance÷burnin) %.\n\nCurrent guess: $p\nJumps = $jumpsize\n")
    flush(stdout)

    @inbounds for i=Base.OneTo(chainsteps)

        ϕ, jumpname, jump = proposaljump(p, jumpsize, rng=rng)
        if strictpriors(ϕ, record_max_age, climate_limits)

            porewaterhistory!(sc, ϕ, k, climate, seawater, ka_dt)

            llCl = loglikelihood(prior.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p) 
            llO = loglikelihood(prior.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)
            llϕ = llCl + llO
            #println("Cl: $llCl, O: $llO") # for troubleshooting.
        else
            llϕ=-Inf
        end

        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            jumpsize = update(jumpsize,jumpname,abs(jump)*scalejump) # update jumpsize
            p = ϕ  # update proposal
            ll = llϕ # Record new log likelihood   
            acceptance[i] = true           
        end

        chains[:,i] .= p.onset, p.dfrz, p.dmlt, p.sea2frz, p.frz2mlt
        lldist[i] = ll

        if iszero(i % chainupdate) # Update progress 
            println("Main Chain --- ", stopwatch(i,chainsteps,clock))
            flush(stdout)
        end

    end
    outnames = (fieldnames(Proposal)...,:ll, :accept)
    outvalues = ((chains[i,:] for i in axes(chains,1))..., lldist, acceptance)
    NamedTuple{outnames}(outvalues), jumpsize
end


