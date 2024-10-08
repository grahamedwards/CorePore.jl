@test Water(1,2).Cl == Water(1.,2.).Cl
@test (Water(1,2)...,) === (1.,2.)
@test mcmurdoshelf().O == -0.33 && mcmurdoshelf().Cl == 19.2657
@test mcmurdosound().O == -1.0 && mcmurdosound().Cl == 19.81655
@test 119.44 < deepbonney().Cl < 119.45 && deepbonney().O ==-25.2

@test isnan(CoreData([1,2], [1,1], [1,1], [1., missing], [1,1]).O.mu[2])
@test CoreData([1,2], [1,1], [1,1], [1., missing], [1,1]).z == [1.0, 2.0]
@test CoreData([1,2], [1,1], [1,1], :Cl).z == [1., 2.]
@test isempty(CoreData([1,2], [1,1], [1,1], :Cl).O.mu)
@test isempty(CoreData([1,2], [1,1], [1,1], :O).Cl.mu)

@test andrill2a().z[1] ≈ 9.67
@test andrill2a().Cl.mu[1] ≈ 23.1843
@test isnan( last( andrill2a().O.mu ) )

@test andrill1b().z[1] ≈ 9.95
@test andrill1b().Cl.mu[1] ≈ 16.845537420070002
@test isnan( last( andrill1b().O.mu ) )

@test length(PorewaterProperty(4).p) == 4
@test PorewaterProperty(2,1).o == PorewaterProperty([1.,1.],[1.,1.]).o

@test length(SedimentColumn(4).Cl.p) == 4
@test SedimentColumn(4, 1., 1.,).rho.o == SedimentColumn(PorewaterProperty(4,1.),PorewaterProperty(4,1.),PorewaterProperty(4,CorePore.density(1.))).rho.o

lr04test = LR04()
@test first(lr04test.t) == 5320.
@test last(lr04test.x) == 3.23

ktest = Constants(k=0.1, dz=5, dt=10, depth=2000)

@test ktest.depth == 2000.
@test ktest.dtdz == 5. * 10.
@test ktest.penultimate_node == ktest.nz-1
@test ktest.nz == 401
@test first(ktest.k1w) ≈ 0.006126451459331848
@test first(ktest.k1cl) ≈ 0.004187765579838535
@test last(ktest.k2cl) ≈ 0.0021567190997326957
@test last(ktest.k2w) ≈ 0.0030695870568939743


proposaltest = Proposal(1,1,1,1,1,1,1,1)
@test update(proposaltest, :onset, 2.) == Proposal(2,1,1,1,1,1,1,1)
@test update(proposaltest, :dfrz, 2.) == Proposal(1,2,1,1,1,1,1,1)
@test update(proposaltest, :dmlt, 2.) == Proposal(1,1,2,1,1,1,1,1)
@test update(proposaltest, :sea2frz, 2.) == Proposal(1,1,1,2,1,1,1,1)
@test update(proposaltest, :frz2mlt, 2.) == Proposal(1,1,1,1,2,1,1,1)
@test update(proposaltest, :flr, 2.) == Proposal(1,1,1,1,1,2,1,1)
@test update(proposaltest, :basalCl, 2.) == Proposal(1,1,1,1,1,1,2,1)
@test update(proposaltest, :basalO, 2.) == Proposal(1,1,1,1,1,1,1,2)

@test (proposaltest...,) == CorePore.fastsplat(proposaltest)


@test CorePore.ProposalPriors(lr04test,ktest).onset === (0.0, 5320.0)
@test CorePore.ProposalPriors(lr04test,ktest, onset=(1.,2.)).onset === (1.,2.)
@test CorePore.ProposalPriors(lr04test,ktest).climatelimits == (2.65, 5.08)
@test CorePore.ProposalPriors(lr04test,ktest, climatelimits=(1.,2.)).climatelimits === (1.,2.)
@test CorePore.ProposalPriors(lr04test,ktest).flr === (0.0, 2000.)
@test CorePore.ProposalPriors(lr04test,ktest, flr=(1,2)).flr === (1.,2.)