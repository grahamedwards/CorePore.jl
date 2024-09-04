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

!!! DEPRECATED !!!
    CorePore.strictpriors(p::Proposal, record_max_age::Number, climatelimits::Tuple{Number,Number}, depth::Number)

Evalute strict constraints on priors that will automatically reject a proposal with..

"""
function strictpriors(p::Proposal, record_max_age::Number, climatelimits::Tuple{Number,Number}, depth::Number)
    
    x=true
    
    x &= 0 < p.onset <= record_max_age
    x &= climatelimits[1] < p.sea2frz < p.frz2mlt < climatelimits[2]
    x &= 0 < p.dfrz <= 0.002 # ≤ 0.4dt/dz to prevent an error in a log-calculation in `CorePore.boundaryconditions`
    x &= 0 < p.dmlt <= 10.
    x &= 0 < p.flr <= depth
    x &= 0 < p.basalCl <= 200
    x &= -56 < p.basalO

    x
end

"""

    checkpriors(p::Proposal, pp::ProposalPriors)

Returns `false` if any proposal value `p` falls beyond the prescribed prior bounds in `pp`, or if proposed subglacial thresholds are non-physical, i.e. melting at lower benthic δ¹⁸O than freezing.
    
Otherwise returns `true`.

"""
function checkpriors(p::Proposal, pp::ProposalPriors)
    
    x=true
    
    x &= pp.onset[1] < p.onset <= pp.onset[2]
    x &= pp.dfrz[1] < p.dfrz <= pp.dfrz[2] 
    x &= pp.dmlt[1] < p.dmlt <= pp.dmlt[2]
    x &= pp.climatelimits[1] < p.sea2frz < p.frz2mlt < pp.climatelimits[2]
    x &= pp.flr[1] < p.flr <= pp.flr[2]
    x &= pp.basalCl[1] < p.basalCl <= pp.basalCl[2]
    x &= pp.basalO[1] < p.basalO < pp.basalO[2]

    x
end



"""

    CorePore.proposaljump(p::Proposal, σ::Proposal; f, rng::AbstractRNG)

Add a random jump to a randomly selected field of `p` with a corresponding normal jumping distribution defined by the corresponding field in `σ`. The possible fields may be specified by providing a Tuple of Symbols `f` (default=`fieldnames(Proposal)`), and a specific RNG seed may be provided.

Note that `:dmlt`, `:dfrz`, and `:basalCl` are drawn from a lognormal jumping distribution (where `σ` is in log-space.)

see also: [`porewatermetropolis`](@ref), [`Proposal`](@ref)

"""
function proposaljump(p::Proposal, j::Proposal; rng::AbstractRNG=Xoshiro(), f::Tuple{Vararg{Symbol}}=fieldnames(Proposal))

    jumpname = rand(rng,f)
    logdist = jumpname ∈ (:dmlt, :dfrz, :basalCl)
    jump = getproperty(j,jumpname) * randn(rng)
    x = getproperty(p,jumpname)
    x = logdist ? log(x) : x
    x += jump
    x = ifelse(logdist, exp(x),x)
    (update(p, jumpname, x) , jumpname , abs(jump) )
end



"""

```julia
porewatermetropolis(p, σ, priors, coredata, cliamte, k; burnin=10, chainsteps=10, climate, k, seawater, onlychloride=true, explore, rng)
```

Executes a Markov chain Monte Carlo (MCMC) routine that explores the parameter space of the variables in [`Proposal`](@ref) instance `p`, constrained by sediment column porewater chemistry records in [`CoreData`](@ref) instance `coredata` and climate record in [`ClimateHistory`](@ref) instance `climate` (e.g. [`LR04`](@ref)), and diffusion-advection parameters in [`Constants`](@ref) instance `k` (I recommend using default values `Constants()`).

`Proposal` instance `σ` describes the (1σ) Gaussian jump size corresponding to each field in `p` (note that fields `:dmlt`, `:dfrz`, and `:basalCl` use a log-normal jump and require a log-space value in `σ`), and `priors` are a [`ProposalPriors`](@ref) instance describing the prior bounds on `p`. 

Kwarg `onlychloride` determines whether the MCMC inverts for only porewater chloridity (default=`true`) or porewater chloridity and δ¹⁸O (`false`).

Additional kwargs include...

kwarg | `DataType` | description | default
:---- | :--------: | :---------- | ------:
burnin | `Int` | number of Markov chain burnin-warm-up steps | `10`
chainsteps | `Int` | number of recorded Markov chains steps | `10`
seawater | [`Water`](@ref) | seawater composition overlying sediment column | `mcmurdosound()`
explore | `Tuple{Vararg{Symbol}}` | parameters explored by the MCMC | `fieldnames(Proposal)`
rng | `AbstractRNG` | optional random number seed | `Random.Xoshiro()`

see also: [`Proposal`](@ref), [`CoreData`](@ref), [`ClimateHistory`](@ref), [`Constants`](@ref), [`Water`](@ref)

"""
function porewatermetropolis(p::Proposal, jumpsigma::Proposal, priors::ProposalPriors, coredata::CoreData, climate::ClimateHistory, k::Constants; burnin::Int=10, chainsteps::Int=10,  seawater::Water=mcmurdosound(), onlychloride::Bool=true, explore::Tuple{Vararg{Symbol}}=fieldnames(Proposal),  rng::AbstractRNG=Random.Xoshiro())

    scalejump=2.4

    pwhfunc = if onlychloride 
        printstyled("Chloride-only mode (onlychloride=true)", bold=true,color=:cyan)
        explore_ = ()
        @inbounds for i ∈ explore
            if i != :basalO
                explore_ = (explore_...,i)
            else
                println(":  MCMC will not explore :basalO (= $(p.basalO))\n\n")
            end
        end
        explore = explore_
        chlorporewaterhistory!
    else
        printstyled("Chloride & δ¹⁸O mode (onlychloride=false)\n\n", bold=true,color=:light_blue)
        porewaterhistory!
    end
    
    flush(stdout)

    ϕ = p # make a new proposal from the original.

    chains = Matrix{Float64}(undef, length(fieldnames(Proposal)), chainsteps)
    lldist = Vector{Float64}(undef, chainsteps)
    acceptance = falses(chainsteps)

    sc = SedimentColumn(k.nz,seawater...)
    ka_dt = CorePore.dt_climatetimestep(climate.t,k.dt)
    
    pwhfunc(sc, ϕ, k, climate, seawater, ka_dt)
    
    llCl, llO = loglikelihood(coredata.z,coredata.Cl.mu,coredata.Cl.sig,k.z,sc.Cl.p), loglikelihood(coredata.z,coredata.O.mu,coredata.O.sig,k.z,sc.O.p)
    ll = llCl + ifelse(onlychloride,0,llO)
    
    clock = time()
    burnupdate, chainupdate = div.((ifelse(iszero(burnin),1,burnin), ifelse(iszero(chainsteps),1,chainsteps)),20,RoundUp)
    println("Beginning sequence...\n  $burnin burn-in iterations \n  $chainsteps recorded iterations\n ------------ \n\n " )
    flush(stdout)

    burninacceptance=0
    
    @inbounds for i=Base.OneTo(burnin)


        ϕ, jumpname, jump = proposaljump(p, jumpsigma, f=explore, rng=rng)
        if checkpriors(ϕ, priors)

            pwhfunc(sc, ϕ, k, climate, seawater, ka_dt)

            llCl, llO = loglikelihood(coredata.z,coredata.Cl.mu,coredata.Cl.sig,k.z,sc.Cl.p), loglikelihood(coredata.z,coredata.O.mu,coredata.O.sig,k.z,sc.O.p)
            llϕ = llCl + ifelse(onlychloride,0,llO)
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
        if checkpriors(ϕ, priors)

            pwhfunc(sc, ϕ, k, climate, seawater, ka_dt)

            llCl, llO = loglikelihood(coredata.z,coredata.Cl.mu,coredata.Cl.sig,k.z,sc.Cl.p), loglikelihood(coredata.z,coredata.O.mu,coredata.O.sig,k.z,sc.O.p)
            llϕ = llCl + ifelse(onlychloride,0,llO)
        else
            llϕ=-Inf
        end
        
        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            jumpsigma = update(jumpsigma, jumpname, jump * scalejump) # update jumpsigma
            p, ll = ϕ, llϕ  # update proposal and log-likelihood
            acceptance[i] = true   
        end
        

        chains[:,i] .= fastsplat(p)
        lldist[i] = ll

        if iszero(i % chainupdate) # Update progress 
            println("Main Chain --- ", stopwatch(i,chainsteps,clock))
            flush(stdout)
        end

    end
    outnames = (fieldnames(Proposal)..., :ll, :accept)
    outvalues = ((chains[i,:] for i in axes(chains,1))..., lldist, acceptance)
    NamedTuple{outnames}(outvalues)
end


