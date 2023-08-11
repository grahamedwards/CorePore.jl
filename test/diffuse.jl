seawater = (cCl =19.2657, d18O = -0.3300)
bctest = (dt=10.,dz=5.,m=2e-4, f=2e-5)
    
    
@test PorewaterDiffusion.boundaryconditions(0., 0., .4, 1., 3., bctest.m, bctest.f, seawater..., bctest.dz, bctest.dt) == (seawater.cCl,seawater.d18O)
    
warmbased = PorewaterDiffusion.boundaryconditions(seawater..., 4., 1., 3., bctest.m, bctest.f, seawater..., bctest.dz, bctest.dt)
    
@test warmbased[1] ≈ 19.246453546453548
@test warmbased[2] ≈ -0.36963036963036966
    
coldbased = PorewaterDiffusion.boundaryconditions(seawater..., 2., 1., 3., bctest.m, bctest.f, seawater..., bctest.dz, bctest.dt)
    
@test coldbased[1] ≈ 19.267626762676265
@test coldbased[2] ≈ -0.33015900795053005


@test PorewaterDiffusion.velocity(2.,4.,.1) ≈ 0.1 
@test PorewaterDiffusion.density(1) ≈ 1.0018



## Test diffusion calculations

k=constants(k=0.1, dz=5, dt=10, depth=2000)
seawater = AND2A()
sc = SedimentColumn(k.nz,seawater.Cl, seawater.O)

sc.Cl.o[1] *= 1.2
sc.O.o[1] *= 1.2
sc.rho.o[1] = density(sc.Cl.o[1])

Cltest = diffusionadvection(sc.Cl.o[2], sc.Cl.o[1], sc.Cl.o[3], k.k1cl[2], k.k2cl[2], PorewaterDiffusion.velocity(sc.rho.o[2],sc.rho.o[1],k.k),k.dt,k.dz) 

d18Otest = diffusionadvection(sc.O.o[2], sc.O.o[1], sc.O.o[3], k.k1w[2], k.k2w[2], PorewaterDiffusion.velocity(sc.rho.o[2],sc.rho.o[1],k.k),k.dt,k.dz) 

diffuseadvectcolumn!(sc,k)

@test Cltest ≈ 19.970844957241017
@test d18Otest ≈ -1.0082012433725027

@test sc.Cl.o[2] == sc.Cl.p[2] ≈ Cltest
@test sc.O.o[2] == sc.O.o[2] ≈ d18Otest