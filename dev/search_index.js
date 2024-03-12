var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PorewaterDiffusion","category":"page"},{"location":"#PorewaterDiffusion","page":"Home","title":"PorewaterDiffusion","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PorewaterDiffusion.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PorewaterDiffusion]","category":"page"},{"location":"#PorewaterDiffusion.ClimateHistory","page":"Home","title":"PorewaterDiffusion.ClimateHistory","text":"ClimateHistory\n\nDataType declared as a shorthand for      NamedTuple{(:t, :x, :n), Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}, Vector{Float64}, Int64}}\n\nThe type generated by LR04, which is used for (perhaps unnecessary) type stability and code readability\n\nsee also: LR04\n\n\n\n\n\n","category":"type"},{"location":"#PorewaterDiffusion.Constants","page":"Home","title":"PorewaterDiffusion.Constants","text":"Constants( ; k, dt, dz, depth )\n\nReturns a Constants struct containing constants and coefficients used in diffusion history calculations. The inputs and their default values are listed in the table below. From these it calculates a few convenient variables–-the product dt * dz (dtdz), a Range of node depths z (in m), the number of nodes nz, and the penultimate_node–-as well as temperature-dependent diffusion coefficients used in diffusionadvection, using the temperature-depth paramterization of Morin+ 2010. These are returned as Vectors of length nz, where each cell corresponds to a depth node in z. The two coefficients are k1 and k2, both of which are follwed with cl or w to denote Cl⁻ or water, respectively: k1cl, k2cl, k1w, and k2w.\n\nfield description default\nk hydraulic conductivity of sediment column 0.1 m / yr\ndt timestep 10 yrs\ndz node spacing 5 m\ndepth sediment column depth 2000 m\n\n\n\n\n\n","category":"type"},{"location":"#PorewaterDiffusion.CoreData","page":"Home","title":"PorewaterDiffusion.CoreData","text":"CoreData(z, mCl, σCl, mO, σO)\n\nReturns an  CoreData struct with sediment core data formatted for porewatermetropolis — (; z, Cl = (mu, sig), O = (mu, sig)). Inputs are Vectors for values of sample depths z (meters below sea floor), measured (mean) chlorinity mCl and 1σ uncertainty σCl, and measured (mean) δ¹⁸O mO and 1σ uncertainties σO. \n\nVectors must all be of the same length. If Cl and δ¹⁸O sampling is not 1:1, use NaN or anything that is not a subtype of Number (e.g. missing, nothing). While NaN is used internally, CoreData does this conversion for you.\n\n\n\nCoreData(z, m, σ, measurement)\n\nSame as above, but for a core with only chlorinity or δ¹⁸O data (the other Vectors are empty to return null values in loglikelihood calculations). Provide sample depths z, measured means m, 1σ uncertainties σ, and the measurement as a symbol, e.g. :Cl or :O. \n\n\n\nPorewaterDiffusion.jl comes loaded with convenient functions to generate data for andrill2a from Tracy+ 2010 (doi:10.1130/G30849.1) and andrill1b from Pompilio+ 2007 (https://digitalcommons.unl.edu/andrillrespub/37). \n\n\n\n\n\n","category":"type"},{"location":"#PorewaterDiffusion.PorewaterProperty","page":"Home","title":"PorewaterDiffusion.PorewaterProperty","text":"PorewaterProperty(n::Int, [, x])\n\nstruct to contain sediment column poperties at each node for a prior timestep o and present timestep p.\n\nConstructor function returns an instance of PorewaterProperty with vectors of length n.  Optionally provide a value x <: Number to fill vectors with (otherwise values are undefined). \n\n\n\n\n\n","category":"type"},{"location":"#PorewaterDiffusion.Proposal","page":"Home","title":"PorewaterDiffusion.Proposal","text":"Proposal\n\nDataType declared as a shorthand for @NamedTuple(onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64) because Graham likes the splatting and name-indexing functionality of NamedTuples and structs don't have them!\n\nsee also: proposal\n\n\n\n\n\n","category":"type"},{"location":"#PorewaterDiffusion.SedimentColumn","page":"Home","title":"PorewaterDiffusion.SedimentColumn","text":"SedimentColumn(n::Int, [, Cl, O])\n\nstruct to contain PorewaterPropertys for the porewater properties of chloridity (Cl), δ¹⁸O (O), and density (rho).\n\nConstructor function returns an instance of SedimentColumn with PorewaterProperty vectors of length n.  Optionally provide values for Cl, O, and rho (otherwise values are undefined). \n\nsee also: PorewaterProperty, density\n\n\n\n\n\n","category":"type"},{"location":"#PorewaterDiffusion.Water","page":"Home","title":"PorewaterDiffusion.Water","text":"Water\n\nDataType declared as a shorthand for NamedTuple{(:Cl, :O), Tuple{Float64,Float64}} because Graham likes the splatting functionality of NamedTuples and structs don't have it!\n\nsee also: water\n\n\n\n\n\n","category":"type"},{"location":"#PorewaterDiffusion.LR04-Tuple{}","page":"Home","title":"PorewaterDiffusion.LR04","text":"LR04()\n\nGenerate a NamedTuple instance containing the Liesiecki & Raymo 2004 benthic stack (data, publication) interpolated for 1 ka timesteps and going forward in time from 5.32 Ma. \n\nfield description\nt time (ka)\nx benthic δ¹⁸O (‰)\nn timesteps in t\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.andrill1b-Tuple{}","page":"Home","title":"PorewaterDiffusion.andrill1b","text":"andrill1b\n\nGenerate a CoreData instance with data from core ANDRILL-1B (from Pompilio+ 2007 (https://digitalcommons.unl.edu/andrillrespub/37)). \n\nsee also: CoreData\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.andrill2a-Tuple{}","page":"Home","title":"PorewaterDiffusion.andrill2a","text":"andrill2a\n\nGenerate a CoreData instance with data from core ANDRILL-2A (from Frank+ 2010, doi:10.1130/G30849.1). \n\nsee also: CoreData\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.boundaryconditions-Tuple{Float64, Float64, Vararg{Any, 9}}","page":"Home","title":"PorewaterDiffusion.boundaryconditions","text":"function boundaryconditions(Cl, d18O, x, sea2freeze,freeze2melt, meltrate, freezerate, Clsw, d18Osw, dz, dt)\n\nCalculates sediment surface boundary condition for δ¹⁸O (d18O) and chloridity (Cl), based on the thermodynamic state described by the current current benthic δ¹⁸O value x and the threshold values corresponding to subglacial freezing sea2freeze and subglacial melting freeze2melt.\n\nFor melting or freezing states, calculates boundary condition from the assumed meltrate, freezerate, timestep dt, length-step dt, and composition of seawater Clsw and d18Osw.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.deepbonney-Tuple{}","page":"Home","title":"PorewaterDiffusion.deepbonney","text":"deepbonney()\n\nGenerate a [`water`](@ref) NamedTuple with >30 m Lake Bonney water compositions. Use for sediment column basal boundary condition. Chloridity (143.333 g/L  ÷ 1.2 kg/L ≈ 119 g/kg) from Angino+ 1963 ([doi:10.1086/626879](https://doi.org/10.1086/626879)) and δ¹⁸O from Matsubaya+ 1979 ([doi:10.1016/0016-7037(79)90042-5](https://doi.org/10.1016/0016-7037(79)90042-5)).\n\nsee also: water.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.density-Tuple{Any}","page":"Home","title":"PorewaterDiffusion.density","text":"density(chlorinity)\n\nCalculates the density of a water parcel with chlorinity in units g/m³ (rather than kg/m³ for convenience with velocity)\n\nsee also: velocity\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.diffuseadvectcolumn!-Tuple{SedimentColumn, Constants, Float64}","page":"Home","title":"PorewaterDiffusion.diffuseadvectcolumn!","text":"diffuseadvectcolumn!(sc, k, flr)\n\nCalculate diffusive and advective transport of chlorinity and isotope-traced water through a sediment column described by properties in k, a NamedTuple generated by the Constants function.\n\nOverwrites the o and p PorewaterProperty fields of sc – a SedimentColumn with pre-existing conditions of Cl, O, and rho in its o fields.\n\nrelies on: velocity, density, diffusionadvection, \n\nsee also: SedimentColumn, Constants \n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.diffusion-NTuple{4, Any}","page":"Home","title":"PorewaterDiffusion.diffusion","text":"diffusion(x,above,below,k)\n\nCalculate the property of a node in a vertical profile given the effect of diffusion, alone. Returns the property given initial values for the node x, the overlying node above, the underlying node below, and (dt/dz²-scaled) diffusion coefficient k.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.diffusionadvection-NTuple{8, Any}","page":"Home","title":"PorewaterDiffusion.diffusionadvection","text":"diffusionadvection(x,above,below,k1,k2,v,dt,dz)\n\nCalculate the property of a node in a vertical profile given the combined effects of diffusion and advection. Returns the property given initial values for the node x, the overlying node above, the underlying node below, (dt/dz-scaled) diffusion coefficients k1 and k2, vertical advection velocity v, timestep dt, and lengthstep dz. Alternatively provide the product of v * dt * dz for a minor speed-up.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.dt_climatetimestep-Tuple{AbstractRange, Number}","page":"Home","title":"PorewaterDiffusion.dt_climatetimestep","text":"dt_climatetimestep(katime,dt)\n\nA helper function to calculate the number of diffusion model timesteps dt (in years) in each timestep of the climate timescale katime (in kiloannum). Returns an integer. \n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.equilibratecolumn!-Tuple{SedimentColumn, @NamedTuple{Cl::Float64, O::Float64}, @NamedTuple{Cl::Float64, O::Float64}, StepRangeLen, Float64}","page":"Home","title":"PorewaterDiffusion.equilibratecolumn!","text":"equilibratecolumn!(sc, seawater, basalwater, z, flr)\n\nCalculate an equilibrium linear profile for all SedimentColumn vectors in sc between a seafloor seawater and basalwater composition, given node depths z and diffusion-dominated column depth flr.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.getproposal-Tuple{@NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}, Symbol}","page":"Home","title":"PorewaterDiffusion.getproposal","text":"getproposal(p::Proposal, s::Symbol)\n\nReturns the value corresponding to the field of Symbol s in Proposal instance p. Use in lieu of getproperty to avoid allocations. \n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.linterp-NTuple{5, Any}","page":"Home","title":"PorewaterDiffusion.linterp","text":"linterp(x, x₁, Δx, y₂, y₁ )\n\nEstimate the value of y corresponding to x given known coordinate (x₁, y₁), value y₂, and Δx= x₂ - x₁. Note that if x==x₁, y= y₁. \n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.loglikelihood-Union{Tuple{T}, Tuple{T, T, T, AbstractRange{Float64}, T}} where T<:Vector{Float64}","page":"Home","title":"PorewaterDiffusion.loglikelihood","text":"loglikelihood( zₒ, μ, σ, zₘ, m )\n\nCalculate the (relative) log-likelihood that model values in m simulated at depth nodes in zₘ were drawn from normally distributed observations with sample depths, mean values, and 1σ uncertainties in corresponding indices of the Vectors zₒ, μₒ, and σ. The simulated value is linearly interpolated between bounding depth nodes in zₘ.\n\nsee also: normll, linterp\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.mcmurdoshelf-Tuple{}","page":"Home","title":"PorewaterDiffusion.mcmurdoshelf","text":"mcmurdoshelf()\n\nGenerate a water NamedTuple with coretop porewater compositions from the McMurdo ice shelf. Pairs with core ANDRILL-1B.\n\nsee also: water, andrill1b\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.mcmurdosound-Tuple{}","page":"Home","title":"PorewaterDiffusion.mcmurdosound","text":"mcmurdosound()\n\nGenerate a water NamedTuple with southern McMurdo Sound water compositions. Pair with core ANDRILL-2A.\n\nsee also: water, andrill2a\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.normll-Tuple{Float64, Float64, Float64}","page":"Home","title":"PorewaterDiffusion.normll","text":"normll(x, μ, σ)\n\nCalculate the (relative) log-likelihood that an observation x was drawn from the normal distribution with mean μ and standard deviation σ. If μ is a NaN (instead of missing for type homogeneity), returns 0.\n\nNote: This function excludes a constant that will not vary for different vlaues of x to speed up calculation in metropolis. To calculate the absolute log-likelihood, take the log of normpdf\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.normpdf-Tuple{Number, Number, Number}","page":"Home","title":"PorewaterDiffusion.normpdf","text":"normpdf(x, μ, σ)\n\nCalculate the probability density of a normal distribution with mean μ and standard deviation σ at the value x.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.porewaterhistory!-Tuple{SedimentColumn, @NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}, Constants, @NamedTuple{t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}, x::Vector{Float64}, n::Int64}, @NamedTuple{Cl::Float64, O::Float64}, Int64}","page":"Home","title":"PorewaterDiffusion.porewaterhistory!","text":"porewaterhistory!(sc::SedimentColumn, p::Proposal, k::Constants, climhist::NamedTuple, seawater::Water, ka_dt::Int)\n\nIn-place version of porewaterhistory, which takes every input as an arg (rather than some defaults as kwargs). It also requires you to provide ka_dt – the number of diffusion timesteps in each thousand-year climate timestep.\n\nsee also: porewaterhistory\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.porewaterhistory-Tuple{@NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}}","page":"Home","title":"PorewaterDiffusion.porewaterhistory","text":"porewaterhistory(p [; k=Constants(), climatehistory=LR04(), seawater=mcmurdosound()])\n\nCalculate the porewater advection-diffusion history of chlorinity and O-isotope-traced water in a sediment column described by properties in k (::Constants) over a given ClimateHistory (LR04 by default) and coretop seawater compositions.\n\nProposal) instance p describes the sensitivity and response of the system to climate fluctuations as recorded in climatehistory.\n\nSee diffuseadvectcolumn! for the underlying diffusion-advection transport calculations.\n\nsee also: porewaterhistory!, Proposal, Constants, LR04, water\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.porewatermetropolis-Tuple{@NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}, @NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}, CoreData}","page":"Home","title":"PorewaterDiffusion.porewatermetropolis","text":"porewatermetropolis...\n\nNot tested, yet...\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.proposal-NTuple{8, Number}","page":"Home","title":"PorewaterDiffusion.proposal","text":"proposal(onset, dfrz, dmlt, sea2frz, frz2mlt, flr, basalCl, basalO)\n\nReturns a NamedTuple (special DataType PorewaterDiffusion.Proposal) with proposal parameters. All inputs must be of type Number (converts to Float64).\n\n\n\nfield description units\nonset onset of model ka\ndfrz freezing rate m/yr\ndmlt melting rate m/yr\nsea2frz Benthic δ¹⁸O threshold for subglacial freezing ‰\nfrz2mlt Benthic δ¹⁸O threshold for subglacial melting ‰\nflr depth of diffusion-dominated porewater column m\nbasalCl chloridity at base of diffusion-dominated column g/kg\nbasalO δ¹⁸O at base of diffusion-dominated column g/kg\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.proposaljump-Tuple{@NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}, @NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}}","page":"Home","title":"PorewaterDiffusion.proposaljump","text":"PorewaterDiffusion.proposaljump(p::Proposal, σ::Proposal; f=proposals, rng::AbstractRNG)\n\nAdd a random jump to a randomly selected field of p with a corresponding normal jumping distribution defined by the corresponding field in σ. The possible fields may be specified by providing a Tuple of Symbols f, and a specific RNG seed may be provided.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.stopwatch-Tuple{Integer, Integer, Number}","page":"Home","title":"PorewaterDiffusion.stopwatch","text":"stopwatch(i, n, t)\n\nConvenience function for [porewatermetropolis] that returns a String reporting the progress at step i for total steps n with start time t (in s since the epoch).\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.strictpriors-Tuple{@NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}, Number, Tuple{Number, Number}, Constants}","page":"Home","title":"PorewaterDiffusion.strictpriors","text":"PorewaterDiffusion.strictpriors(p::Proposal, record_max_age::Number, climatelimits::Tuple{Number,Number}, k::Constants)\n\nEvalute strict constraints on priors that will automatically reject a proposal with...\n\nOnset date beyond the climate record timespan (in ka) (record_max_age)\nNonphysical subglacial thresholds – melting at lower benthic δ¹⁸O than freezing or values exceeding the record extrema (climatelimits).\nFreezing rate is non-zero and ≤ 0.4dt/dz (to prevent an error in a log-calculation in PorewaterDiffusion.boundaryconditions).\nAnnual melting rate is non-zero and must be no more than that observed at the Thwaites grounding line (<10 m/yr, Davis+ 2023).\nDiffusive porewater column (p.flr) is between 0 and 2 km depth.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.update-Tuple{@NamedTuple{onset::Float64, dfrz::Float64, dmlt::Float64, sea2frz::Float64, frz2mlt::Float64, flr::Float64, basalCl::Float64, basalO::Float64}, Symbol, Number}","page":"Home","title":"PorewaterDiffusion.update","text":"update(p::Proposal, f::Symbol, x::Number)\n\nUpdate field f of Proposal instance p with value x.\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.velocity-Tuple{Any, Any, Any}","page":"Home","title":"PorewaterDiffusion.velocity","text":"velocity(x, above, k)\n\nCalculate the velocity (m/yr) at a node with density x, given the density of the node above, and the hydraulic conductivity k (m/yr).\n\nsee also: density\n\n\n\n\n\n","category":"method"},{"location":"#PorewaterDiffusion.water-Tuple{Number, Number}","page":"Home","title":"PorewaterDiffusion.water","text":"water(Cl, O)\n\nReturns a NamedTuple with values of chlorinity Cl and δ¹⁸O O(necessarily Float64s). \n\nsee also: mcmurdoshelf, mcmurdosound, PorewaterDiffusion.Water\n\n\n\n\n\n","category":"method"}]
}
