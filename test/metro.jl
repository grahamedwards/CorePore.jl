
# strictpriors
spargs = (1000., (1.,4.))
@test PorewaterDiffusion.strictpriors(Proposal(500,1,1,2,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(1100,1,1,2,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(500,-1,1,2,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(500,1,-1,2,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(500,1,1,0,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(500,1,1,2,5), spargs...)


# proposaljump
@test PorewaterDiffusion.proposaljump(Proposal(ones(5)...), Proposal(ones(5)...), rng=StableRNG(2580))[1].dfrz ≈ 1.090668193354693

# stopwatch
@test "0% |■■□□□□□□□□| 100%  ||  step: 27 / 100  ||  time: 0.0 m" == PorewaterDiffusion.stopwatch(27,100,time())