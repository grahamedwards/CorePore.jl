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
Calculate the (relative) log-likelihood that an observation `x` was drawn from the normal distribution with mean `μ` and standard deviation `σ`. If `μ` is a NaN (instead of missing for type homogeneity), returns `0`.

Note: This function excludes a constant that will not vary for different vlaues of `x` to speed up calculation in `metropolis`. To calculate the absolute log-likelihood, take the `log` of `normpdf`
"""
normll(μ::Float64,σ::Float64,x::Float64) = ifelse(isnan(μ),zero(x), -(x-μ)*(x-μ) / (2*σ*σ))





"""

    linterp(x, x₁, Δx, y₂, y₁ )

Estimate the value of y corresponding to `x` given known coordinate (`x₁`, `y₁`), value `y₂`, and `Δx`= x₂ - x₁. Note that if `x`==`x₁`, y= `y₁`. 

"""
linterp(x,x1,dx,y2,y1) = y1 + (x-x1) * (y2-y1)/dx






"""
    loglikelihood( zₒ, μ, σ, zₘ, m ) 

Calculate the (relative) log-likelihood that model values in `m` simulated at depth nodes in `zₘ` were drawn from normally distributed observations with sample depths, mean values, and 1σ uncertainties in corresponding indices of the Vectors `zₒ`, `μₒ`, and `σ`. The simulated value is linearly interpolated between bounding depth nodes in `zₘ`.

see also: [`normll`](@ref), [`linterp`](@ref)

"""
function loglikelihood(zo::T, muo::T, sigo::T, zm::AbstractRange{Float64}, m::T) where T<:Vector{Float64}
    ll=zero(eltype(m))
    dzm = step(zm)

    @inbounds @simd for i = eachindex(zo)
        zoi = zo[i]
        mi = searchsortedfirst(zm, zoi) # find first index ≥ zoi
        mi_ = ifelse(zoi == zm[mi], mi, mi-1) # if zm[mi] ≠ zoi, mi_ is the preceding index, otherwise mi = mi_ and xmi = m[mi]
        xmi = linterp(zoi, zm[mi_], dzm, m[mi], m[mi_])
        ll += normll(muo[i], sigo[i], xmi)
    end
    ll
end


"""

    means(x; std::Integer)

Calculate the means of each field in NamedTuple `x`. Optionally provide an integer number of standard deviations to calculate, returned as the form `key` = (m=μ, s=σ). Use with the returned chains of [`porewatermetropolis`](@ref).

"""
function means(x::NamedTuple; std::Int=0)
    k = keys(x)
    if std > 0
    NamedTuple{k}( (m=Statistics.mean(x[i]), s=std*Statistics.std(x[i])) for i in k) 
    else 
    NamedTuple{k}(Statistics.mean(x[i]) for i in k) 
    end 
end




"""

    medians(x; ci::Float64)

Calculate the medians of each field in NamedTuple `x`. Optionally provide a `ci` ∈ [0,1] of standard deviations to calculate, returned as the form `key` = (m=median, l=lower, u=upper). Use with the returned chains of [`porewatermetropolis`](@ref).

"""
function medians(x; ci::Float64=0.)
    @assert 0 ≤ ci ≤ 1
    k = keys(x)
    if ci>0.
    upper, lower = (0.5 + 0.5ci), 0.5ci
    NamedTuple{k}( (m=Statistics.median(x[i]), l=Statistics.quantile(x[i],lower), u=Statistics.quantile(x[i],upper)) for i in k)
    else 
    NamedTuple{k}(Statistics.median(x[i]) for i in k) 
    end 
end