function cCl = cCl_diffusion(cCl_,cCl_above,cCl_below,Diff_Cl,Diff_Cl_,v,dt,dz)
% CCL_DIFFUSION calculates cCl in a vertical profile
% 
% Inputs: cCl for the current cell in the previous timestep (cCl_),
% cCl of the overlying cell in the previous timestep (cCl_a), 
% cCl of the underlying cell in the previous timestep (cCl_b), 
% diffusion coefficient of Cl at corresponding depth/temperature (Diff_Cl), 
% diffusion coefficient of Cl of the overlying water parcel (Diff_Cl_), 
% vertical velocity (advection) of the porewaters in that parcel, 
% time step (dt), depth step (dz)
% 
% Outputs: chloride concentration of current cell
% 

coeff_Cl = Diff_Cl.*dt./dz./dz; % coefficient used in FD calculations
coeff2_Cl = (Diff_Cl - Diff_Cl_) * dt/dz; % coefficient used in FD calculations

cCl=cCl_ + coeff_Cl*(cCl_above - 2*cCl_ + cCl_below) + coeff2_Cl*(cCl_above - cCl_) - (v*dt*dz)*(cCl_ - cCl_above); % Cl concetration of cell

end