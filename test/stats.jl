

@test CorePore.linterp(2,1,2,4,8) == CorePore.linterp(2,1,2,8,4) == 6.
@test CorePore.linterp(1,1,2,8,4) == 4.


@test CorePore.normll(0.,1.,1.) ≈ -0.5
@test CorePore.normll(NaN,1.,1.) == 0.
@test CorePore.normpdf(0,1,0.5) ≈ 0.3520653267642995

lltest = (; zo=[0, 1.8, 2.2], muo=[.2, .42, .6], sigo = [.08, .12, .14] , zm=0:4., m= [.2,.4, .5, .4, .3] )

@test -0.50 < CorePore.loglikelihood(lltest...) < -0.49

@test iszero(CorePore.loglikelihood(zeros(0), lltest.muo, lltest.sigo, lltest.zm, lltest.m))

#(CorePore.normll(lltest.muo[1],lltest.sigo[1],lltest.m[1]), CorePore.normll(lltest.muo[2],lltest.sigo[2],CorePore.linterp(lltest.zo[2], lltest.zm[3],1, lltest.m[4], lltest.m[3])), CorePore.normll(lltest.muo[3],lltest.sigo[3],CorePore.linterp(lltest.zo[3], lltest.zm[4],1, lltest.m[5],lltest.m[4])))