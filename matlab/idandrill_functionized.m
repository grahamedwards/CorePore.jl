%% 1D ionic diffusion in a sedimentary column
% same as above, but this time with functions


addpath('functions/')

tic

% choose core to examine
% 1 = AND-1B, 2 = AND-2A
core=2;

% choose parameters based on which core you're examining
if core==1 
    seawater.cCl=19.2657; % g/kg, average value of seawater (Andrill MIS data)
    seawater.d18O=-0.33; % per mil, Andrill MIS data

    f_rt=2e-5; % m/yr, freeze-on rate
    m_rt=2e-4; % m/yr, melt rate
end

if core==2
    seawater.cCl=19.81655; % g/kg, average value of seawater (Andrill SMS data)
    seawater.d18O=-1; % per mil, Andrill SMS data

    f_rt=4e-5; % m/yr, freeze-on rate
    m_rt=1e-2; % m/yr, melt rate
end



% load in benthic stack data
% use for setting boundaries for seawater/melting/freezing
Benthic_interp = load ('../data/LR04-interpolated-1ka.csv');

% flip matrix so that the first element is 5320 ka
Benthic=flipud(Benthic_interp); 
% no need to have Benthic_rev and Benthic, just say Benthic=flipup(Benthic)

% set limits separating melting/freezing/ocean conditions
freeze2melt=4.2; % upper limit separating meltwater from freezing
ocean2freeze=3.5; % lower limit separating seawater from freezing

% based on von Neumann stability condition (1/2 >= D dt/(dx^2)):
dt=10; % yrs
dz=5; % m

Depth = 2000; % depth of the domain in meters
nz = Depth/dz + 1; % number of vertical nodes


% number of nodes for time
nStart=length(Benthic(:,1));

% timestep
num=1000/dt; % years

k = 0.1; % m/yr, hydraulic conductivity of sedimentary column 

% starting concentration is concentration of Cl in seawater (McM Sound)
C_Cl = zeros(nz,2); % column vectors containing values of concentration
C_d18O = zeros(nz,2); % column vectors containing values of concentration
% column 2 = new, column 1 = old


% when considering diffusion to be temperature dependent:
depth_vec=0:dz:Depth;
temperature = (0.0767.*depth_vec)-2.4+273.15; % Morin et al. 2010
% Diff_Cl=exp(3.8817+(-2.2854e+03./T)); % Cl
Diff_Cl=exp(3.8817+(-2.2854e+03./temperature)); % Cl, m^2/yr
Diff_d18O=exp(4.2049+(-2.2699e+03./temperature)); % d18O, m^2/yr

% coeff and coeff2 are built into the diffusion functions
coeff_Cl = Diff_Cl.*dt./dz./dz; % coefficient used in FD calculations
coeff_d18O = Diff_d18O.*dt./dz./dz; % coefficient used in FD calculations

coeff2_Cl = zeros(1,length(coeff_Cl));
coeff2_d18O = zeros(1,length(coeff_d18O));

for i=2:length(Diff_Cl)
    coeff2_Cl(i) = (Diff_Cl(i) - Diff_Cl(i-1)) * dt/dz;
    coeff2_d18O(i) = (Diff_d18O(i) - Diff_d18O(i-1)) * dt/dz;
end


% initial conditions - assumed that initially all seawater
C_Cl(:,:)=seawater.cCl; 
C_d18O(:,:)=seawater.d18O;

% rho=(seawater.cCl * 0.0018 ) + 1 ; % initial porewater density


% The main time loop solves the 1D equation: dC/dt = Diff*d2C/dz2 using the
% simplest explicit Euler Finite Difference approximation



for n=1:nStart-1 % timespan of model run

% pick out boundary conditions 
mof=meltoceanfreeze(Benthic(n,2),ocean2freeze,freeze2melt);


% run time loop
for s = 1:num 
    
% set upper boundary conditions
[cClo,d18Oo]=chooseboundaries(mof,C_Cl(1,1),C_d18O(1,1),dz,dt,m_rt,f_rt,seawater);
C_Cl(1,1) = cClo; % setting the upper boundary condition
C_Cl(1,2) = cClo; % setting the upper boundary condition
C_d18O(1,1) = d18Oo; % setting the upper boundary condition
C_d18O(1,2) = d18Oo; % setting the upper boundary condition

  

% account for density driven vertical flow of brine in sediments
    rho=((C_Cl(:,1).*0.0018)+1)*1000; % fluid density
    if rho(i-1)<rho(i)
        v=0;
    else
        v=k.*(rho(i-1)-rho(i))./rho(i); % vertical velocity
    end    


% run diffusion model
for i = 2:(nz-1)
    
%     [rho,v]=rhov(rho,C_Cl(i,1),k);

    C_Cl(i,2) = cCl_diffusion(C_Cl(i,1),C_Cl(i-1,1),C_Cl(i+1,1),coeff_Cl(i),coeff2_Cl(i),v,dt,dz);
    C_d18O(i,2) = d18O_diffusion(C_d18O(i,1),C_d18O(i-1,1),C_d18O(i+1,1),coeff_d18O(i),coeff2_d18O(i),v,dt,dz);

end
C_Cl(end,2)=C_Cl(end-1,2); % set bottom value to penultimate value
C_d18O(end,2)=C_d18O(end-1,2); % set bottom value to penultimate value

C_Cl(:,1) = C_Cl(:,2); % putting the new values into the 'old' vector
C_d18O(:,1) = C_d18O(:,2); % putting the new values into the 'old' vector


end % time loop ends 



% n


end % end run



toc
