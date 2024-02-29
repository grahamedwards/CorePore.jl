@test PorewaterDiffusion.dt_climatetimestep(100:-1:0,10) == 100

ch, k, sw, p = LR04(), Constants(), mcmurdosound(), Proposal(5320., 4e-5,1e-2,3.5, 4.2)
sc = SedimentColumn(k.nz,sw...)
porewaterhistory!(sc,p, k, ch,sw, sw, PorewaterDiffusion.dt_climatetimestep(ch.t,k.dt))

@test 108 < sc.Cl.p[end] < 109
@test 19 < sc.Cl.p[2] < 20

@test -5.3 < sc.O.p[end] < -5.2
@test -2.6 < sc.O.p[2] < -2.5

pwh2 = porewaterhistory(p,k=k,climatehistory=ch,seawater=sw,basalwater=sw)

@test 108 < pwh2.Cl[end] < 109
@test 19 < pwh2.Cl[2] < 20

@test -5.3 < pwh2.d18O[end] < -5.2
@test -2.6 < pwh2.d18O[2] < -2.5