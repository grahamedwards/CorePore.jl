% Very helpful link: 
% https://www.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
% 
% tips:
%   assert(test, 'message')      ... e.g. assert(x > 0, 'Must be positive')
%   testresult = runtests('test/tests.m'); % runs tests
%   table(testresult) % returns a nice table of test results.

% Unit tests
% includes:
%   chooseboundaries

%% chooseboundaries
% the %% above means this will be reported in the table as a block. Most excellent! 

cbtol = 1e-4; % set some reasonable accuracy tolerance (I

seawater.cCl = 19.2657; % build a struct with seawater compositions.
seawater.d18O = -0.3300; 

cClo = 0.; %set as zero so we can test the seawater condition (mof=0)
d18Oo = 0.;
dt = 10.;
dz = 5.;
meltingrate=2e-4;
freezingrate=2e-5;

% run the function in each unique scenario:
[cClsw,d18Osw] = chooseboundaries(0,cClo,d18Oo, dz, dt, meltingrate, freezingrate, seawater);
[cClf,d18Of] = chooseboundaries(1,cClsw,d18Osw, dz, dt, meltingrate, freezingrate, seawater);
[cClm,d18Om] = chooseboundaries(-1,cClsw,d18Osw, dz, dt, meltingrate, freezingrate, seawater);


assert(cClsw == seawater.cCl, 'seawater cClo failed')
assert(d18Osw == seawater.d18O, 'seawater d18Oo failed')
assert(abs(19.2676-cClf)<cbtol,'freezing condition cClo failed') % give a meaningful message.
assert(abs(-0.3302-d18Of)<cbtol,'freezing condition d18Oo failed') 

assert(abs(-0.3696-d18Om)<cbtol,'melting condition d18Oo failed') 
assert(abs(19.2465-cClm)<cbtol,'melting condition cClo failed')

%% rhov
rhovtest_rho = 2.;
rhovtest_cCl = 20.;
rhovtest_k = 0.1;

[rho,v] = rhov(rhovtest_rho,rhovtest_cCl,rhovtest_k)
assert(v>0,'velocity not positive when rho_<rho')
assert(abs(rho-1.0359)<1e-4,'vertical velocity calculation failed')

[rho,v] = rhov(rho,rhovtest_cCl,rhovtest_k)
assert(v==0, 'velocity non-zero when rho_ >= rho')

%% nextfunction