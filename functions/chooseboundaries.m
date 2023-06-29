function [cClo,d18Oo] = chooseboundaries(mof,cClo,d18Oo, dz, dt, m_rt, f_rt, seawater)
% CHOOSEBOUNDARIES calculate the boundary conditions
%   Returns the surface boundary conditions of Cl- concentration cClo and
%   delta-18O (d18Oo) for three conditions:
%   seawater -> mof = 0
%   melting -> mof = -1
%   freezing -> mof = 1
%
% Takes as inputs coretop [Cl-] (cClo) and delta-18O (d18Oo), lengthscale
% (dz), timestep (dt), melting rate (m_rt), freezing rate (f_rt), and
% seawater composition (struct containing fields cCl and d18O).
%   
% see also MELTOCEANFREZ % this could be a simple function that returns a
% -1, 0, or 1 for an input of the current benthic value and the threshold
% variables.

if mof==0 % ocean conditions
    cClo = seawater.cCl;
    d18Oo = seawater.d18O;

elseif mof==1 % freezing conditions
    dz_ = dz*0.4;
    frtdt = f_rt*dt;

    cClo = cClo * dz_ / (dz_ - frtdt); 
    d18Oo = d18Oo + 1.59 * log(1 - frtdt / dz_); % simplified from eqn 2 of Toyota et al. (2017)
   

elseif mof==-1 % melting conditions

    dz_ = dz*0.4;
    mrtdt = m_rt*dt;
    
    cClo = cClo * dz_ / (dz_ + mrtdt); 
    d18Oo =  (d18Oo *  dz_ - 40 * mrtdt) / (dz_+mrtdt);
end

end