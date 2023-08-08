seawater = (cCl =19.2657, d18O = -0.3300)
bctest = (dt=10.,dz=5.,m=2e-4, f=2e-5)
    
    
@test PorewaterDiffusion.boundaryconditions(0., 0., .4, 1., 3., bctest.m, bctest.f, seawater..., bctest.dz, bctest.dt) == (seawater.cCl,seawater.d18O)
    
warmbased = PorewaterDiffusion.boundaryconditions(seawater..., 4., 1., 3., bctest.m, bctest.f, seawater..., bctest.dz, bctest.dt)
    
@test warmbased[1] ≈ 19.246453546453548
@test warmbased[2] ≈ -0.0006593406593406596
    
coldbased = PorewaterDiffusion.boundaryconditions(seawater..., 2., 1., 3., bctest.m, bctest.f, seawater..., bctest.dz, bctest.dt)
    
@test coldbased[1] ≈ 19.267626762676265
@test coldbased[2] ≈ -0.33015900795053005


@test PorewaterDiffusion.velocity(2.,4.,.1) ≈ 0.1 
@test PorewaterDiffusion.density(1) ≈ 1.0018