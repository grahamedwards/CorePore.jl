

@test PorewaterDiffusion.linterp(2,1,2,4,8) == PorewaterDiffusion.linterp(2,1,2,8,4) == 6.
@test PorewaterDiffusion.linterp(1,1,2,8,4) == 4.


@test PorewaterDiffusion.normll(0.,1.,1.) ≈ -0.5
@test PorewaterDiffusion.normll(NaN,1.,1.) == 0.
@test PorewaterDiffusion.normpdf(0,1,0.5) ≈ 0.3520653267642995

lltest = (; zo=[0, 1.8, 2.2], muo=[.2, .42, .6], sigo = [.08, .12, .14] , zm=0:4, m= [.2,.4, .5, .4, .3] )

@test -0.72 <PorewaterDiffusion.loglikelihood(lltest...) < -0.71 

PorewaterDiffusion.loglikelihood(Vector{Float64}(undef,0), lltest.muo, lltest.sigo, lltest.zm, lltest.m)
#(PorewaterDiffusion.normll(lltest.muo[1],lltest.sigo[1],lltest.m[1]), PorewaterDiffusion.normll(lltest.muo[2],lltest.sigo[2],PorewaterDiffusion.linterp(lltest.zo[2], lltest.zm[3],1, lltest.m[4], lltest.m[3])), PorewaterDiffusion.normll(lltest.muo[3],lltest.sigo[3],PorewaterDiffusion.linterp(lltest.zo[3], lltest.zm[4],1, lltest.m[5],lltest.m[4])))