sw = mcmurdoshelf()
bctest = (dt=10.,dz=5.,m=2e-4, f=2e-5)
    
    
@test PorewaterDiffusion.boundaryconditions(0., 0., .4, 1., 3., bctest.m, bctest.f, sw..., bctest.dz, bctest.dt) == (sw...,1.03467826)
    
warmbased = PorewaterDiffusion.boundaryconditions(sw..., 4., 1., 3., bctest.m, bctest.f, sw..., bctest.dz, bctest.dt)
    
@test warmbased[1] ≈ 19.246453546453548
@test warmbased[2] ≈ -0.36963036963036966
    
coldbased = PorewaterDiffusion.boundaryconditions(sw..., 2., 1., 3., bctest.m, bctest.f, sw..., bctest.dz, bctest.dt)
    
@test coldbased[1] ≈ 19.267626762676265
@test coldbased[2] ≈ -0.33015900795053005


@test PorewaterDiffusion.velocity(2.,4.,.1) ≈ 0.1 
@test PorewaterDiffusion.density(1) ≈ 1.0018



## Test diffusion calculations

k=constants(k=0.1, dz=5, dt=10, depth=2000)
sw = mcmurdosound()
sc = SedimentColumn(k.nz,sw...)

sc.Cl.o[1]= sc.Cl.o[end] *= 1.2
sc.O.o[1] = sc.O.o[end] *= 1.2
sc.rho.o[1] = density(sc.Cl.o[1])

Cltest = diffusionadvection(sc.Cl.o[2], sc.Cl.o[1], sc.Cl.o[3], k.k1cl[2], k.k2cl[2], PorewaterDiffusion.velocity(sc.rho.o[2],sc.rho.o[1],k.k),k.dt,k.dz) 

d18Otest = diffusionadvection(sc.O.o[2], sc.O.o[1], sc.O.o[3], k.k1w[2], k.k2w[2], PorewaterDiffusion.velocity(sc.rho.o[2],sc.rho.o[1],k.k),k.dt,k.dz) 

@test Cltest ≈ 19.970844957241017
@test d18Otest ≈ -1.0082012433725027

diffuseadvectcolumn!(sc,k)

@test sc.Cl.o == sc.Cl.p
@test sc.O.o == sc.O.p
@test sc.rho.o == sc.rho.p

@test sc.Cl.o[2] ≈ Cltest
@test sc.Cl.o[end] ≈ 20.166286404580042
@test sc.O.o[2] ≈ d18Otest
@test sc.O.o[end] ≈ -1.025290758205530