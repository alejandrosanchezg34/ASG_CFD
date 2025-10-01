% ----------------------------------------------------------------------- %
% Test case: Implementing HGS to solve a nozzle flow
% ----------------------------------------------------------------------- %
% Author: Sanchez Garcia, Alejandro
% Master final thesis 2025
% ESEIAAT UPC
% ----------------------------------------------------------------------- %

close all
clear
clc

addpath("..\");
addpath("..\HGS_main\");

% ----------------------------------------------------------------------- %

%% 1 - Input data

% Select Riemann solver:
Riemann_solver = "HLLC"; % Either "HLLC" or "exact"

% Configuration:
% - Subsonic inflow: fix inlet density and pressure, but let inlet velocity
% float.
% - Supersonic outflow: let all variables float (density, velocity and
% pressure).
BCs = "subsonic_in_free_out";

% Numerical data:
N = 50; % Number of grid nodes
sigma = 0.9; % Stability criteria (CFL)
t_end = 25E-2; % [s] - Desired final time
w = 0; % [] - Constant parameter for MUSCL-Hancock's slope vectors (can go from -1 to 1)

% Geometric data:
L = 3; % [m] - Domain length
At = 0.01; % [m2] - Nozzle throat area

% Gas properties:
gas_model = "frozen_flow"; % Select gas model: "constant", "frozen_flow" or "shifting_flow".
species = {'H2';'O2';'H2O';'H';'O';'OH'}; % Species involved in the gas
rOF = 8; % O2/H2 ratio (stoichiometric) at combustion chamber
Tr = 298; % [K] - Temperature of reactants
p_0 = 80; % [bar] - Combustion chamber pressure
mH2 = 1/(1+rOF); % H2 mass fraction
mO2 = 1 - mH2; % O2 mass fraction
MH2 = HGSsingle('H2','Mm')/1000; % [kg/mol] - H2 molecular mass
MO2 = HGSsingle('O2','Mm')/1000; % [kg/mol] - O2 molecular mass
nH2 = mH2/MH2; % [mols] - H2
nO2 = mO2/MO2; % [mols] - O2
nr = [nH2;nO2;0;0;0;0]; % [mols] - Composition of reactants

% Get combustion chamber temperature and products composition:
[T_0,~,n_0] = HGStp(species,nr,'T',Tr,p_0); % [K], [mols]

% Get the gas properties of the combustion chamber products:
[Rg_0,gamma_0] = HGSprop(species,n_0,T_0,p_0,'Rg','gamma'); % [kJ/(kgK)], [] - Gas constant and adiabatic coefficient
Rg_0 = Rg_0*1E3; % [J/(kgK)] - Gas constant
cv_0 = getcv(Rg_0,gamma_0); % [J/(kgK)] - Gas specific heat at constant volume

% Flow properties:
p_0 = p_0*1E5; % [Pa] - Stagnation pressure
rho_0 = p_0/(Rg_0*T_0); % [kg/m3] - Stagnation density
a_0 = sqrt(gamma_0*Rg_0*T_0); % [m/s] - Stagnation sound speed
p_out = NaN; % [Pa] - Exit pressure (not needed)

%% 2 - 1D mesh generation

% Generating mesh:
[x,dx] = generate_mesh(L,N,"Godunov"); % [m], [m]

% Generating nozzle cross-section area:
% A_f = @(x) (x <= L/2)*(At*(1 + 0.2223*(3*x/L - 1.5)^2)) + ...
%     (x > L/2)*(At*(1 + 2.2*(3*x/L - 1.5)^2)); % [m2] - Cross-section area analytical function
% dAdx_f = @(x) (x <= L/2)*(At*0.2223*2*3/L*(3*x/L - 1.5)) + ...
%     (x > L/2)*(At*2.2*2*3/L*(3*x/L - 1.5)); % Cross-section area derivative analytical function
A_f = @(x) At*(1 + 2.2*(3*x/L - 1.5)^2); % [m2] - Cross-section area analytical function
dAdx_f = @(x) At*2.2*2*3/L*(3*x/L - 1.5); % Cross-section area derivative analytical function

A = arrayfun(A_f,x); % [m2] - Cross-section area distribution vector

%% 3 - Initial conditions

% Initial conditions (at t = 0):
rho(1,1:N) = rho_0; % [kg/m3] - Density
T(1,1:N) = T_0; % [K] - Temperature
u(1,1:N) = 0.10*a_0; % [m/s] - Velocity
e(1,1:N) = cv_0*T(1,:); % [m2/s2] - Internal specific energy
p(1,1:N) = p_0; % [Pa] - Pressure

%% 4 - Solver

% Preparing flow-field variables structure:
vars.rho = rho; % [kg/m3] - Density
vars.u = u; % [m/s] - Velocity
vars.e = e; % [m2/s2] - Specific energy
vars.T = T; % [K] - Temperature
vars.p = p; % [Pa] - Pressure

% Solving:
tic
[vars,gas,tvec] = SOLVER_MUSCL_Hancock(vars,A_f,dAdx_f,x,dx,t_end,sigma,w,rho_0,NaN,p_0,p_out,BCs,Riemann_solver, ...
    gas_model,Rg_0,gamma_0,species,n_0);
solver_time = toc;
fprintf("Computation time: %.4f s \n",solver_time);

% Decoding flow-field variables from solver structure:
rho = vars.rho; % [kg/m3] - Density
u = vars.u; % [m/s] - Velocity
e = vars.e; % [m2/s2] - Specific energy
T = vars.T; % [K] - Temperature
p = vars.p; % [Pa] - Pressure

% Decoding gas properties from solver structure:
Rg = gas.Rg; % [J/(kgK)] - Gas constant
gamma = gas.gamma; % [] - Gas adiabatic coefficient
cv = getcv(Rg,gamma); % [J/(kgK)] - Gas specific heat at constant volume coefficient

%% 5 - Postprocess

% Compute other interesting data (from numerical solution):
a = sqrt(gamma.*Rg.*T); % [m/s] - Sound speed
M = u./a; % [] - Mach number
mdot = rho.*A.*u; % [kg/s] - Mass flow

% Compute last time step composition (if shifting flow):
if gas_model == "shifting_flow"
    n_end = zeros(length(species),N);
    for i = 1:N
        [~,n_end(:,i)] = HGSeq(species,n_0,T(end,i),p(end,i)*1E-5); % [mols] - Equilibrium composition
        n_end(:,i) = n_end(:,i)./sum(n_end(:,i))*100; % Convert to percentage
    end
end

%% 6 - Plotting variable distributions

figure

% Temperature:
subplot(2,2,1)
hold on
plot(x,T(end,:),"-o","DisplayName",gas_model)
hold off
xlabel("x [m]")
ylabel("T [K]")
grid on
legend("Location","best")
title("Temperature")

% Pressure:
subplot(2,2,2)
hold on
plot(x,p(end,:),"-o","DisplayName",gas_model)
hold off
xlabel("x [m]")
ylabel("p [Pa]")
grid on
legend("Location","best")
title("Pressure")

% Velocity:
subplot(2,2,3)
hold on
plot(x,u(end,:),"-o","DisplayName",gas_model)
hold off
xlabel("x [m]")
ylabel("u [m/s]")
grid on
legend("Location","best")
title("Velocity")

% Mach number:
subplot(2,2,4)
hold on
plot(x,M(end,:),"-o","DisplayName",gas_model)
hold off
xlabel("x [m]")
ylabel("Mach")
grid on
legend("Location","best")
title("Mach number")

%% 7 - Plotting gas properties and composition

% Composition:
figure
hold on
for j = 1:length(species)
    if gas_model == "shifting_flow"
        plot(x,n_end(j,:),"-o","DisplayName",species{j,1})
    else
        plot(x,(n_0(j)/sum(n_0))*ones(1,N)*100,"-o","DisplayName",species{j,1})
    end
end
hold off
xlabel("x [m]")
ylabel("Molar composition [%]")
grid on
legend("Location","best")
title("Gas composition")

figure

% Gas constant:
subplot(2,2,[1 2])
hold on
plot(x,Rg(end,:),"-o","DisplayName",gas_model)
hold off
xlabel("x [m]")
ylabel("R_g [J/(kgK)]")
grid on
legend("Location","best")
title("Gas constant")

% Gas adiabatic coefficient:
subplot(2,2,3)
hold on
plot(x,gamma(end,:),"-o","DisplayName",gas_model)
hold off
xlabel("x [m]")
ylabel("\gamma")
grid on
legend("Location","best")
title("Gas adiabatic coefficient")

% Gas adiabatic coefficient:
subplot(2,2,4)
hold on
plot(x,cv(end,:),"-o","DisplayName",gas_model)
hold off
xlabel("x [m]")
ylabel("c_v [J/(kgK)]")
grid on
legend("Location","best")
title("Gas specific heat")

% ----------------------------------------------------------------------- %