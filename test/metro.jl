
# strictpriors
spargs = (1000., (1.,4.), Constants())
@test PorewaterDiffusion.strictpriors(proposal(500,.001,1,2,3,1,1,1), spargs...) # passing test
@test !PorewaterDiffusion.strictpriors(proposal(1100,.001,1,2,3,1,1,1), spargs...) # onset too old
@test !PorewaterDiffusion.strictpriors(proposal(0,.001,1,2,3,1,1,1), spargs...) # onset too young
@test !PorewaterDiffusion.strictpriors(proposal(500,.001,1,0,3,1,1,1), spargs...) # sea2frz too low
@test !PorewaterDiffusion.strictpriors(proposal(500,.001,1,2,5,1,1,1), spargs...) # frz2mlt too high
@test !PorewaterDiffusion.strictpriors(proposal(500,.001,1,3,2,1,1,1), spargs...) # sez2frz < frz2mlt 
@test !PorewaterDiffusion.strictpriors(proposal(500,.01,1,2,3,1,1,1), spargs...) # p.dfrz too high.
@test !PorewaterDiffusion.strictpriors(proposal(500,.001,20.,2,3,1,1,1), spargs...) # p.dfrz too high.
@test !PorewaterDiffusion.strictpriors(proposal(500,.001,1,2,3,2001,1,1), spargs...) # p.flr too deep
@test !PorewaterDiffusion.strictpriors(proposal(500,.001,1,2,3,0,1,1), spargs...) # p.flr too deep

# proposaljump
@test PorewaterDiffusion.proposaljump(proposal(ones(8)...), proposal(ones(8)...), rng=StableRNG(2580))[1].dfrz ≈ 1.090668193354693 # if lognormal -> 1.0949056480096304

@test PorewaterDiffusion.proposaljump(proposal(ones(8)...), proposal(ones(8)...), rng=StableRNG(10))[1].onset ≈ 0.23064139237391934

# stopwatch
@test "0% |■■□□□□□□□□| 100%  ||  step: 27 / 100  ||  time: 0.0 m" == PorewaterDiffusion.stopwatch(27,100,time())