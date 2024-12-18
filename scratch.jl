using CorePore, Serialization;

cd(@__DIR__)
filename="240611"

println(filename)

p = Proposal(5320., 1e-4,1e-3,3.5, 4.2, 1000, mcmurdoshelf()...)
jumpsize = Proposal(20., .1e-4, .1e-3, .1, .1,10, 10, 1.)

x = andrill1b()

prp = (:dfrz, :dmlt, :sea2frz, :frz2mlt, :flr, :basalCl)
println(prp)

chains = porewatermetropolis(p, jumpsize, x; burnin=20, chainsteps=20, k=Constants(), seawater=mcmurdoshelf(), climate=LR04(),explore=prp, onlychloride=true)

serialize(filename*"run",chains)

println("Serialized chains successfully saved to $filename.cereal")

exit()



using CairoMakie

prior = andrill2a()

x = CorePore.sieveresults(chains, k=Constants(), climate=LR04(), seawater = mcmurdosound(), sieve=1, start=1, stop=0)
