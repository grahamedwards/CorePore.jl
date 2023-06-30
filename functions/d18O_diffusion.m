function d18O = d18O_diffusion(d18O_,d18O_above,d18O_below,Diff_d18O,Diff_d18O_,v,dt,dz)
% CCL_DIFFUSION calculates cCl in a vertical profile
% 
% Inputs: d18O for the current cell in the previous timestep (d18O_),
% d18O of the overlying cell in the previous timestep (d18O_above), 
% d18O of the underlying cell in the previous timestep (cCl_below), 
% diffusion coefficient of d18O at corresponding depth/temperature (Diff_d18O), 
% diffusion coefficient of d18O of the overlying water parcel (Diff_d18O_), 
% vertical velocity (advection) of the porewaters in that parcel, 
% time step (dt), depth step (dz)
% 
% Outputs: d18O of current cell
% 

coeff_d18O = Diff_d18O.*dt./dz./dz; % coefficient used in FD calculations
coeff2_d18O = (Diff_d18O - Diff_d18O_) * dt/dz; % coefficient used in FD calculations

d18O=d18O_ + coeff_d18O*(d18O_above - 2*d18O_ + d18O_below) + coeff2_d18O*(d18O_above - d18O_) - (v*dt*dz)*(d18O_ - d18O_above); % d18O of cell

end