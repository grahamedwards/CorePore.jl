# CorePore

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://grahamedwards.github.io/CorePore.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://grahamedwards.github.io/CorePore.jl/dev/)
[![Build Status](https://github.com/grahamedwards/CorePore.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/grahamedwards/CorePore.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/grahamedwards/CorePore.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/grahamedwards/CorePore.jl)

Diffusive-advective evolution of porewater chemistry in a marine sediment column with changing seafloor chemistry linked to climate and ice sheet grounding line behavior.

Based on diffusion-advection model framework of [Neuhaus & Tulaczyk 2023](https://doi.org/10.1017/aog.2023.28), translated into [Julia](https://julialang.org/) from the original MATLAB for speed 🚀

The model is wrapped in a Markov chain Monte Carlo algorithm. Given a climate record, the MCMC inverts porewater chemistry data for posterior chains of the onset of glaciation at the marine site, the climate sensitivity of the system, and the basal boundary conditions of the sediment column.

## Installation

To install `CorePore.jl` on your own computer, just type `]` into the Julia REPL to enter the built-in package manager and then type:
```julia
add https://github.com/grahamedwards/CorePore.jl
```
 and hit enter.

After installing, just type `using CorePore` to use. 

## Usage

The primary functionality of this package is hosted in the function `porewatermetropolis`. Please read the documentation for its full functionality and the underlying forward models. The example below will give you a quick start:

```julia
using CorePore

p = Proposal(5320., 1e-4,1e-3,3.5, 4.2, 1000, deepbonney()...)
jumpsize = Proposal(20., .1e-4, .1e-3, .1, .1,10, 10, 1.)

chains = porewatermetropolis(p, jumpsize, andrill2a(); burnin=1000, chainsteps=1000, k=Constants(), seawater=mcmurdosound(), climate=LR04(), onlychloride=true)
```