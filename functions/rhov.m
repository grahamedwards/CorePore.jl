function [rho,v] = rhov(rho_,cCl,k)
%% RHOV calculate vertical velocity 
%   Inputs: the returned density of the overlying water parcel (rho_), the chlorinity of
%   the current cell, and the hydraulic conductivity k.
% 
%   Outputs: the density of the current cell and its vertical velocity (v)
%

rho=(cCl * 0.0018 ) + 1 ; % fluid density ( /1000, since the 1000 divides out in the v calculation
if rho_ < rho
    v=0;
else
    v= k * (rho_-rho)/rho; % vertical velocity
end

end