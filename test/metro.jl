
# strictpriors

pp = ProposalPriors(LR04(),Constants())
spargs = (1000., (1.,4.), Constants().depth)
@test CorePore.checkpriors(Proposal(500,.001,1,3,4,1,1,1), pp) # passing test
@test !CorePore.checkpriors(Proposal(50000,.001,1,3,4,1,1,1), pp) # onset too old
@test !CorePore.checkpriors(Proposal(-2,.001,1,3,4,1,1,1), pp) # onset too young
@test !CorePore.checkpriors(Proposal(500,.001,1,0,4,1,1,1), pp) # sea2frz too low
@test !CorePore.checkpriors(Proposal(500,.001,1,3,8,1,1,1), pp) # frz2mlt too high
@test !CorePore.checkpriors(Proposal(500,.001,1,3,2,1,1,1), pp) # sez2frz < frz2mlt 
@test !CorePore.checkpriors(Proposal(500,.01,1,3,4,1,1,1), pp) # p.dfrz too high.
@test !CorePore.checkpriors(Proposal(500,.001,20.,3,4,1,1,1), pp) # p.dfrz too high.
@test !CorePore.checkpriors(Proposal(500,.001,1,3,4,2001,1,1), pp) # p.flr too deep
@test !CorePore.checkpriors(Proposal(500,.001,1,3,4,0,1,1), pp) # p.flr too deep
@test !CorePore.checkpriors(Proposal(500,.001,1,3,4,1,-1,1), pp) # p.basalCl too low
@test !CorePore.checkpriors(Proposal(500,.001,1,3,4,1,201,1), pp) # p.basalCl too high
@test !CorePore.checkpriors(Proposal(500,.001,1,3,4,1,1,-60), pp) # p.basal too low

# proposaljump
@test CorePore.proposaljump(Proposal(ones(8)...), Proposal(ones(8)...), rng=StableRNG(2580))[1].dfrz ≈ 1.0949056480096304 # if linear-sapce -> 1.090668193354693

@test CorePore.proposaljump(Proposal(ones(8)...), Proposal(ones(8)...), rng=StableRNG(10))[1].onset ≈ 0.23064139237391934

# stopwatch
@test "0% |■■□□□□□□□□| 100%  ||  step: 27 / 100  ||  time: 0.0 m" == CorePore.stopwatch(27,100,time())



p, jumpsize = Proposal(5320., 1e-4,1e-3, 3.5, 4.2, 1000., deepbonney()...,), Proposal(10., .1, .1, .1, .1,10, 1,1)

climate, k = LR04(),Constants()
pp = ProposalPriors(climate, k)

@silence x = porewatermetropolis(p, jumpsize, pp, andrill2a(), climate, k; burnin=20, chainsteps=20, seawater=mcmurdosound(), onlychloride=false, rng=StableRNG(2560))

@test first(x.ll) < last(x.ll)
@test first(x.basalCl) > last(x.basalCl)
@test first(x.basalO) < last(x.basalO)
@test first(x.onset) > last(x.onset)

@silence x = porewatermetropolis(p, jumpsize, pp, andrill2a(), climate, k; burnin=20, chainsteps=20, seawater=mcmurdosound(), onlychloride=true, rng=StableRNG(2560))

@test first(x.ll) < last(x.ll)
@test last(x.basalO) == p.basalO
