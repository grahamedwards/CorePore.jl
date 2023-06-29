%% 1D ionic diffusion in a sedimentary column
% loosely based on equations and approach from Adkins and Schrag, 2003,
% Chloride, and d18O
% MIS Core, SMS Core
% benthic stack
% cleaner version
tic

% choose core to examine
% 1 = AND-1B, 2 = AND-2A
core=2;


% load in benthic stack data
% use for setting boundaries for seawater/melting/freezing
Benthic=readmatrix('LR04stack.csv','range','A6:B2120');

% flip matrix so that the first element is 5320 ka
Benthic_rev=flipud(Benthic);


% set boundary condition handle
upper=4.2; % upper limit separating meltwater from freezing
lower=3.5; % lower limit separating seawater from freezing
bc=zeros(length(Benthic_rev(:,2)),1);
bc(Benthic_rev(:,2) < lower)=1; % ocean conditions
bc(Benthic_rev(:,2) < upper & Benthic_rev(:,2) >= lower)=2; % freezing conditions
bc(Benthic_rev(:,2) >= upper)=3; % melting conditions


% number of nodes for time
nStart=length(Benthic_rev(:,1));

% choose parameters based on which core you're examining
if core==1 
    Cl_McM=19.2657; % g/kg, average value of seawater (Andrill MIS data)
    d18O_McM=-0.33; % per mil, Andrill MIS data

    f_rt=2e-5; % m/yr, freeze-on rate
    m_rt=2e-4; % m/yr, melt rate
end

if core==2
    Cl_McM=19.81655; % g/kg, average value of seawater (Andrill SMS data)
    d18O_McM=-1; % per mil, Andrill SMS data

    f_rt=4e-5; % m/yr, freeze-on rate
    m_rt=1e-2; % m/yr, melt rate
end




% based on von Neumann stability condition (1/2 >= D dt/(dx^2)):
dt=10; % yrs
dz=5; % m

Depth = 2000; % depth of the domain in meters
nz = Depth/dz + 1; % number of vertical nodes



k = 0.1; % m/yr, hydraulic conductivity of sedimentary column 

% starting concentration is concentration of Cl in seawater (McM Sound)
C_Cl = zeros(nz,2); % column vectors containing values of concentration
C_d18O = zeros(nz,2); % column vectors containing values of concentration
% column 2 = new, column 1 = old

% need to figure out rate of Cl incorporation to the surface
% 19.81655 g/kg = Cl of seawater above AND-2A core, data repository for Frank et al 2010
% Cl_McM=19.81655; % g/kg
% % figure out how much Cl excluded by basal freezing
% m=1e-3; % m/yr, freeze-on rate, number we've calculated
% % assume 1 m^2 area, McM seawater freezing, only freshwater is freezing
% % Cl_addition = Cl_McM * m * (1 m^2) * rho_water
% Cl_addition=Cl_McM*m*1000; % g/yr
% % need to add this to topmost element (bin) in model
% % volume of water in bin = 1m2 * dz * porosity
% % Cl_input = Cl_addition mass / mass_water = Cl_addition / (rho_water * (1m2) * dz * porosity)
% Cl_input=Cl_addition/(1000*dz*0.4); % g/kg /yr






% when considering diffusion to be temperature dependent:
depth_vec=0:dz:Depth;
temperature = (0.0767.*depth_vec)-2.4+273.15; % Morin et al. 2010
% Diff_Cl=exp(3.8817+(-2.2854e+03./T)); % Cl
Diff_Cl=exp(3.8817+(-2.2854e+03./temperature)); % Cl, m^2/yr
Diff_d18O=exp(4.2049+(-2.2699e+03./temperature)); % d18O, m^2/yr

coeff_Cl = Diff_Cl.*dt./dz./dz; % coefficient used in FD calculations
coeff_d18O = Diff_d18O.*dt./dz./dz; % coefficient used in FD calculations

coeff2_Cl = zeros(1,length(coeff_Cl));
coeff2_d18O = zeros(1,length(coeff_d18O));

for i=2:length(Diff_Cl)
    coeff2_Cl(i) = (Diff_Cl(i) - Diff_Cl(i-1)) * dt/dz;
    coeff2_d18O(i) = (Diff_d18O(i) - Diff_d18O(i-1)) * dt/dz;
end



N_Cl=zeros(nz,1);
N_d18O=zeros(nz,1);

rho=zeros(nz,1);
v=zeros(nz,1);

% initial conditions - assumed that initially all seawater
C_Cl(:,:)=Cl_McM; 
C_d18O(:,:)=d18O_McM; 


% The main time loop solves the 1D equation: dC/dt = Diff*d2C/dz2 using the
% simplest explicit Euler Finite Difference approximation



for n=1:nStart-1 % timespan of model run


    % years spent in timestep
    time=(Benthic_rev(n,1)-Benthic_rev(n+1,1))*1000; % years 


% run time loop
for s = 1:time/dt 
%% This whole block is now a function ~ `chooseboundaries`
% pick out boundary conditions and set upper boundary conditions
    
if bc(n)==1 % ocean conditions
    C_Cl(1,1) = Cl_McM; % setting the upper boundary condition
    C_Cl(1,2) = Cl_McM; % setting the upper boundary condition
    C_d18O(1,1) = d18O_McM; % setting the upper boundary condition
    C_d18O(1,2) = d18O_McM; % setting the upper boundary condition

else 
if bc(n)==2 % freezing conditions

    % figure out how much Cl or d18O excluded by basal freezing
    % f_rt=5e-5; % m/yr, freeze-on rate
    % assume 1 m^2 area, only freshwater is freezing
    % new concentration = mass Cl/(mass water originally in bin - mass water freezing)
    % Cl_new = C_Cl(1,1)*Vbin/(Vbin-Vloss)
    % Cl_new = C_Cl(1,1)*dz*porosity/((dz*porosity)-(f_rt*dt))
    C_Cl(1,1) = C_Cl(1,1)*dz*0.4/((dz*0.4)-(f_rt*dt)); % setting the upper boundary condition
    C_Cl(1,2) = C_Cl(1,2)*dz*0.4/((dz*0.4)-(f_rt*dt)); % setting the upper boundary condition

    % eqn 2 from Toyota et al. (2017)
    % C_new-C_initial=1.59*log(V_ice/V_water_initial)
    C_d18O(1,1)=C_d18O(1,1)+ 1.59 * log((dz*.4-f_rt*dt)/(dz*.4)); % setting the upper boundary condition
    C_d18O(1,2)=C_d18O(1,2)+ 1.59 * log((dz*.4-f_rt*dt)/(dz*.4)); % setting the upper boundary condition

else
if bc(n)==3 % melting conditions

    % figure out how much H2O added by basal melting
    % m_rt=1e-3; % m/yr, melt rate
    % assume 1 m^2 area
    % Cl_concentration = C_Cl(1,1)*Vw_bin/(Vw_bin + Vw_add)
    % Vw_bin = dz * (1m2) * porosity, Vw_add = m * (1m2) * time
    % Cl=C_Cl*dz*0.4/((dz*0.4)+(m_rt*dt));
    
    C_Cl(1,1) = C_Cl(1,1)*dz*0.4/((dz*0.4)+(m_rt*dt)); % setting the upper boundary condition
    C_Cl(1,2) = C_Cl(1,2)*dz*0.4/((dz*0.4)+(m_rt*dt)); % setting the upper boundary condition
    % setting the upper boundary condition
    C_d18O(1,1) = C_d18O(1,1)*((dz*0.4)/((dz*0.4)+(m_rt*dt))) + (-40)*((m_rt*dt)/((dz*0.4)+(m_rt*dt)));
    C_d18O(1,2) = C_d18O(1,2)*((dz*0.4)/((dz*0.4)+(m_rt*dt))) + (-40)*((m_rt*dt)/((dz*0.4)+(m_rt*dt)));

end
end
end
    
%% The next lines rho= ... end would make a great function. 
% To speed it up, rather than calculate a vector of rho, you could make a
% function that goes into the for-loop below. This function takes a value rho and a value
% from C_Cl, and returns both a new rho (to be reused in the next loop) and
% a value for v. (I think the calculation is faster for the whole vector,
% but reallocating the new vector is gonna ding you every time)

% account for density driven vertical flow of brine in sediments
    rho=((C_Cl(:,1).*0.0018)+1)*1000; % fluid density
    if rho(i-1)<rho(i)
        v=0;
    else
        v=k.*(rho(i-1)-rho(i))./rho(i); % vertical velocity
    end

%% Turn each diffusion calculation into a function.
% run diffusion model
for i = 2:(nz-1)
    
    % diffusion is depth-dependent
    C_Cl(i,2) = C_Cl(i,1) + coeff_Cl(i)*(C_Cl(i-1,1)-2*C_Cl(i,1)+C_Cl(i+1,1)) + coeff2_Cl(i)*(C_Cl(i-1,1)-C_Cl(i,1)) - (v*dt*dz)*(C_Cl(i,1)-C_Cl(i-1,1));
    C_d18O(i,2) = C_d18O(i,1) + coeff_d18O(i)*(C_d18O(i-1,1)-2*C_d18O(i,1)+C_d18O(i+1,1)) + coeff2_d18O(i)*(C_d18O(i-1,1)-C_d18O(i,1)) - (v*dt*dz)*(C_d18O(i,1)-C_d18O(i-1,1));
end
C_Cl(end,2)=C_Cl(end-1,2); % set bottom value to penultimate value
C_d18O(end,2)=C_d18O(end-1,2); % set bottom value to penultimate value

C_Cl(:,1) = C_Cl(:,2); % putting the new values into the 'old' vector
C_d18O(:,1) = C_d18O(:,2); % putting the new values into the 'old' vector


end % time loop ends 



% n

N_Cl(:)=C_Cl(:,2);
N_d18O(:)=C_d18O(:,2);
end % end run

toc