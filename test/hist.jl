@test PorewaterDiffusion.dt_climatetimestep(100:-1:0,10) == 100

ch, k, sw = LR04(), Constants(), mcmurdosound()
p = Proposal(5320., 4e-5,1e-2,3.5, 4.2, 1990., sw...)
sc = SedimentColumn(k.nz,sw...)
porewaterhistory!(sc,p, k, ch,sw, PorewaterDiffusion.dt_climatetimestep(ch.t,k.dt))

@test sc.Cl.p[end-1] == sw.Cl
@test 19 < sc.Cl.p[2] < 20

@test sc.O.p[end-1] == sw.O
@test -2.6 < sc.O.p[2] < -2.5

sc = SedimentColumn(k.nz,sw...)
PorewaterDiffusion.chlorporewaterhistory!(sc,p, k, ch,sw, PorewaterDiffusion.dt_climatetimestep(ch.t,k.dt))
@test 19 < sc.Cl.p[2] < 20
@test sc.O.p[2] == sw.O

pwh2 = porewaterhistory(p,k=k,climatehistory=ch,seawater=sw)

@test  pwh2.Cl[end-1] == sw.Cl
@test 19 < pwh2.Cl[2] < 20

@test pwh2.d18O[end-1] == sw.O
@test -2.6 < pwh2.d18O[2] < -2.5