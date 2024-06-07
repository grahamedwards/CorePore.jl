
# strictpriors
spargs = (1000., (1.,4.), Constants().depth)
@test PorewaterDiffusion.strictpriors(Proposal(500,.001,1,2,3,1,1,1), spargs...) # passing test
@test !PorewaterDiffusion.strictpriors(Proposal(1100,.001,1,2,3,1,1,1), spargs...) # onset too old
@test !PorewaterDiffusion.strictpriors(Proposal(0,.001,1,2,3,1,1,1), spargs...) # onset too young
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,1,0,3,1,1,1), spargs...) # sea2frz too low
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,1,2,5,1,1,1), spargs...) # frz2mlt too high
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,1,3,2,1,1,1), spargs...) # sez2frz < frz2mlt 
@test !PorewaterDiffusion.strictpriors(Proposal(500,.01,1,2,3,1,1,1), spargs...) # p.dfrz too high.
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,20.,2,3,1,1,1), spargs...) # p.dfrz too high.
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,1,2,3,2001,1,1), spargs...) # p.flr too deep
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,1,2,3,0,1,1), spargs...) # p.flr too deep
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,1,2,3,1,-1,1), spargs...) # p.basalCl too low
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,1,2,3,1,201,1), spargs...) # p.basalCl too high
@test !PorewaterDiffusion.strictpriors(Proposal(500,.001,1,2,3,1,1,-60), spargs...) # p.basal too low

# proposaljump
@test PorewaterDiffusion.proposaljump(Proposal(ones(8)...), Proposal(ones(8)...), rng=StableRNG(2580))[1].dfrz ≈ 1.0949056480096304 # if linear-sapce -> 1.090668193354693

@test PorewaterDiffusion.proposaljump(Proposal(ones(8)...), Proposal(ones(8)...), rng=StableRNG(10))[1].onset ≈ 0.23064139237391934

# stopwatch
@test "0% |■■□□□□□□□□| 100%  ||  step: 27 / 100  ||  time: 0.0 m" == PorewaterDiffusion.stopwatch(27,100,time())



p, jumpsize = Proposal(5320., 1e-4,1e-3, 3.5, 4.2, 1000., deepbonney()...,), Proposal(10., .1, .1, .1, .1,10, 1,1)



@silence x = porewatermetropolis(p, jumpsize, andrill2a(); burnin=20, chainsteps=20, k=Constants(), seawater=mcmurdosound(), climate=LR04(), onlychloride=false, rng=StableRNG(2560))

@test first(x.ll) < last(x.ll)
@test first(x.basalCl) > last(x.basalCl)
@test first(x.basalO) < last(x.basalO)
@test first(x.onset) > last(x.onset)

@silence x = porewatermetropolis(p, jumpsize, andrill2a(); burnin=20, chainsteps=20, k=Constants(), seawater=mcmurdosound(), climate=LR04(), onlychloride=true, rng=StableRNG(2560))

@test first(x.ll) < last(x.ll)
@test last(x.basalO) == p.basalO
