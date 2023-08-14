@test PorewaterDiffusion.dt_climatetimestep(100:-1:0,10) == 100

ch = LR04()
k = constants()
sw = AND2A()
kadt = PorewaterDiffusion.dt_climatetimestep(ch.t,k.dt)
p = Proposal(5320., 4e-5,1e-2,3.5, 4.2)

sc = SedimentColumn(k.nz,sw...)
pwh = porewaterhistory!(sc,p, k, ch,sw, kadt)

@test 108 < pwh.Cl[end] < 109
@test 19 < pwh.Cl[2] < 20

@test -5.3 < pwh.d18O[end] < -5.2
@test -2.6 < pwh.d18O[2] < -2.5