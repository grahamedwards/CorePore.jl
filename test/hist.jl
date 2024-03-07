@test PorewaterDiffusion.dt_climatetimestep(100:-1:0,10) == 100

ch, k, sw = LR04(), Constants(), mcmurdosound()
p = proposal(5320., 4e-5,1e-2,3.5, 4.2, 2000, sw...)
sc = SedimentColumn(k.nz,sw...)
porewaterhistory!(sc,p, k, ch,sw, PorewaterDiffusion.dt_climatetimestep(ch.t,k.dt))

@test 48 < sc.Cl.p[end-1] < 52
@test 19 < sc.Cl.p[2] < 20

@test -2.0 < sc.O.p[end-1] < -1.9
@test -2.6 < sc.O.p[2] < -2.5

pwh2 = porewaterhistory(p,k=k,climatehistory=ch,seawater=sw)

@test 48 < pwh2.Cl[end-1] < 52
@test 19 < pwh2.Cl[2] < 20

@test -2 < pwh2.d18O[end-1] < -1.9
@test -2.6 < pwh2.d18O[2] < -2.5