
# strictpriors
spargs = (1000., (1.,4.), Constants())
@test PorewaterDiffusion.strictpriors(Proposal(500,.1,1,2,3), spargs...) # passing test
@test !PorewaterDiffusion.strictpriors(Proposal(1100,.1,1,2,3), spargs...) # onset too old
@test !PorewaterDiffusion.strictpriors(Proposal(0,.1,1,2,3), spargs...) # onset too young
@test !PorewaterDiffusion.strictpriors(Proposal(500,.1,1,0,3), spargs...) # sea2frz too low
@test !PorewaterDiffusion.strictpriors(Proposal(500,.1,1,2,5), spargs...) # frz2mlt too high
@test !PorewaterDiffusion.strictpriors(Proposal(500,.1,1,3,2), spargs...) # sez2frz < frz2mlt 
@test !PorewaterDiffusion.strictpriors(Proposal(500,.3,1,2,3), spargs...) # p.dfrz too high.
@test !PorewaterDiffusion.strictpriors(Proposal(500,.1,20.,2,3), spargs...) # p.dfrz too high.


# proposaljump
@test PorewaterDiffusion.proposaljump(Proposal(ones(5)...), Proposal(ones(5)...), rng=StableRNG(2580))[1].dfrz ≈ 1.0949056480096304

@test PorewaterDiffusion.proposaljump(Proposal(ones(5)...), Proposal(ones(5)...), rng=StableRNG(10))[1].onset ≈ 0.23064139237391934

# stopwatch
@test "0% |■■□□□□□□□□| 100%  ||  step: 27 / 100  ||  time: 0.0 m" == PorewaterDiffusion.stopwatch(27,100,time())