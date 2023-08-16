
# strictprios
spargs = (1000., (1.,4.))
@test PorewaterDiffusion.strictpriors(Proposal(500,1,1,2,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(1100,1,1,2,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(500,-1,1,2,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(500,1,-1,2,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(500,1,1,0,3), spargs...)
@test !PorewaterDiffusion.strictpriors(Proposal(500,1,1,2,5), spargs...)


# proposaljump
pjtest = PorewaterDiffusion.proposaljump(Proposal(ones(5)...), Proposal(ones(5)...), rng=StableRNG(2580))
pjtest[1].dfrz ≈ 1.090668193354693


# stopwatch
"0% |■■□□□□□□□□| 100%  ||  total: 0.0 m  ||  step: 27 / 100\n" == PorewaterDiffusion.stopwatch(27,100,time())