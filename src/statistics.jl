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

    linterp(x, y₂, y₁, x₁, Δx)

Estimate the value of y corresponding to `x` given known coordinate (`x₁`, `y₁`), value `y₂`, and `Δx`= x₂ - x₁. Note that if `x`==`x₁`, y= `y₁`. 

"""
linterp(x,x1,dx,y2,y1) = y1 + (x-x1) * (y2-y1)/dx