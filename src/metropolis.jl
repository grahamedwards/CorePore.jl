

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
Not tested, not exported... yet...
"""
function porewatermetropolis(p::Proposal, jumpsize::Proposal, prior::NamedTuple; burnin::Int=0, chainsteps::Int=100, k::Constants=Constants(), seawater::Seawater=mcmurdosound(), benthic::CoreData=LR04(), scalejump=2.9, rng::AbstractRNG=Random.Xoshiro())

    record_max_age = first(benthic.t)
    benthic_limits = extrema(benthic.x)

    ϕ = p # make a new proposal from the original.

    chains = Matrix{Float64}(undef, length(fieldnames(Proposal)), chainsteps)
    lldist = Vector{Float64}(undef, chainsteps)
    acceptance = falses(chainsteps)

    sc = SedimentColumn(k.nz,seawater...)
    ka_dt = PorewaterDiffusion.dt_climatetimestep(benthic.t,k.dt)
    
    porewaterhistory!(sc, ϕ, k, benthic, seawater, ka_dt)

    ll= loglikelihood(prior.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p) + loglikelihood(prior.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)

    clock = time()
    println("Beginning sequence...\n  $burnin burn-in iterations \n  $chainsteps recorded iterations\n ------------ \n\n " )
    flush(stdout)
    @inbounds for i=Base.OneTo(burnin)

        ϕ, jumpname,jump = proposaljump(p, jumpsize,rng=rng)
        if strictpriors(ϕ, record_max_age, benthic_limits)

            porewaterhistory!(sc, ϕ, k, benthic, seawater, ka_dt)

            llϕ = loglikelihood(prior.z,prior.Cl.mu,prior.Cl.sig,k.z,sc.Cl.p) + loglikelihood(prior.z,prior.O.mu,prior.O.sig,k.z,sc.O.p)
        else
            llϕ=-Inf
        end

        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            jumpsize = update(jumpsize,jumpname,abs(jump)*scalejump) # update jumpsize
            p = ϕ  # update proposal
            ll = llϕ # Record new log likelihood              
        end

        if iszero(i % 500) # Update progress every 500 steps
            println("Burn-In --- ", stopwatch(i,burnin,clock))
            flush(stdout)
        end
    end

    println(burnin," burn-in steps complete. Current guess: ",p,", ℓ = $ll" )
    flush(stdout)

    @inbounds for i=Base.OneTo(chainsteps)

        ϕ, jumpname,jump = proposaljump(p, jumpsize)
        if strictpriors(ϕ, record_max_age, benthic_limits)

            porewaterhistory!(sc, ϕ, k, benthic, seawater, ka_dt)

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

        if iszero(i % 500) # Update progress every 500 steps
            println("Main Chain --- ", stopwatch(i,chainsteps,clock))
            flush(stdout)
        end

    end
    outnames = (fieldnames(Proposal)...,:ll, :accept)
    outvalues = ((chains[i,:] for i in axes(chains,1))..., lldist, acceptance)
    NamedTuple{outnames}(outvalues)
end


