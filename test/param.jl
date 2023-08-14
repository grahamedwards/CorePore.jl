@test seawater(1,2).Cl == seawater(1.,2.).Cl
@test AND1B().O == -0.33 && AND1B().Cl == 19.2657
@test AND2A().O == -1.0 && AND2A().Cl == 19.81655

@test length(PorewaterProperty(4).p) == 4
@test PorewaterProperty(2,1).o == PorewaterProperty([1.,1.],[1.,1.]).o

@test length(SedimentColumn(4).Cl.p) == 4
@test SedimentColumn(4, 1., 1.,).rho.o == SedimentColumn(PorewaterProperty(4,1.),PorewaterProperty(4,1.),PorewaterProperty(4,PorewaterDiffusion.density(1.))).rho.o

lr04test = LR04()
@test first(lr04test.t) == 5320.
@test last(lr04test.x) == 3.23

ktest = constants(k=0.1, dz=5, dt=10, depth=2000)

@test ktest.depth == 2000.
@test ktest.dtdz == 5. * 10.
@test ktest.penultimate_node == ktest.nz-1
@test ktest.nz == 401
@test first(ktest.k1w) ≈ 0.006126451459331848
@test first(ktest.k1cl) ≈ 0.004187765579838535
@test last(ktest.k2cl) ≈ 0.0021567190997326957
@test last(ktest.k2w) ≈ 0.0030695870568939743