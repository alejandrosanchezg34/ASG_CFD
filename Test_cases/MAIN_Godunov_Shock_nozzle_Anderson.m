% ----------------------------------------------------------------------- %
% Test case: Capturing a shock in Anderson's nozzle
% ----------------------------------------------------------------------- %
% Author: Sanchez Garcia, Alejandro
% Master final thesis 2025
% ESEIAAT UPC
% ----------------------------------------------------------------------- %

close all
clear
clc

addpath("..\");

% ----------------------------------------------------------------------- %

%% 1 - Input data

% Select solver/scheme to use:
solver = "Godunov"; % Either "Godunov" or "MUSCL_Hancock"

% Select Riemann solver:
Riemann_solver = "HLLC"; % Either "HLLC" or "exact"

% Configuration:
% - Subsonic inflow: fix inlet density and pressure, but let inlet velocity
% float.
% - Subsonic outflow: fix outlet pressure, and let outlet density and
% velocity float.
BCs = "subsonic_in_fixed_out";

% Numerical data:
N = 60; % Number of grid nodes
sigma = 0.5; % Stability criteria (CFL)
t_end = 10E-2; % [s] - Desired final time
w = 0; % [] - Constant parameter for MUSCL-Hancock's slope vectors (can go from -1 to 1)

% Geometric data:
L = 3; % [m] - Domain length
At = 0.01; % [m2] - Nozzle throat area

% Gas properties:
gas_model = "constant"; % Constant gas model: properties and composition do not depend on T
species = NaN; % Gas species (NOT NEEDED for the constant gas model)
n_0 = NaN; % Gas composition (NOT NEEDED for the constant gas model)
Rg_0 = 287; % [J/(kgK)] - Gas constant
gamma_0 = 1.4; % [] - Gas adiabatic coefficient
cv_0 = getcv(Rg_0,gamma_0); % [J/(kgK)] - Gas specific heat at constant volume

% Flow properties:
T_0 = 3500; % [K] - Stagnation temperature
p_0 = 1e7; % [Pa] - Stagnation pressure
rho_0 = p_0/(Rg_0*T_0); % [kg/m3] - Stagnation density
a_0 = sqrt(gamma_0*Rg_0*T_0); % [m/s] - Stagnation sound speed
p_out = 0.6784*p_0; % [Pa] - Exit pressure (inner shock wave)

%% 2 - 1D mesh generation

% Generating mesh:
[x,dx] = generate_mesh(L,N,"Godunov"); % [m], [m]

% Generating nozzle cross-section area:
A_f = @(x) At*(1 + 2.2*(3*x/L - 1.5)^2); % [m2] - Cross-section area analytical function
dAdx_f = @(x) At*2.2*2*3/L*(3*x/L - 1.5); % Cross-section area derivative analytical function

A = arrayfun(A_f,x); % [m2] - Cross-section area distribution vector

%% 3 - Initial conditions

% Initializing variable matrices:
rho = zeros(1,N);
u = zeros(1,N);
p = zeros(1,N);
T = zeros(1,N);
e = zeros(1,N);

% Initial conditions (at t = 0):
for i = 1:N
    if (x(i) <= L/6)
        rho(1,i) = rho_0; % [kg/m3] - Density
        T(1,i) = T_0; % [K] - Temperature

    elseif (L/6 < x(i) && x(i) <= L/2)
        rho(1,i) = rho_0*(1.0 - 0.366*(3*x(i)/L - 0.5)); % [kg/m3] - Density
        T(1,i) = T_0*(1.0 - 0.167*(3*x(i)/L - 0.5)); % [K] - Temperature

    elseif (L/2 < x(i) && x(i) <= 0.7*L)
        rho(1,i) = rho_0*(0.634 - 0.702*(3*x(i)/L - 1.5)); % [kg/m3] - Density
        T(1,i) = T_0*(0.833 - 0.4908*(3*x(i)/L - 1.5)); % [K] - Temperature

    else
        rho(1,i) = rho_0*(0.5892 + 0.10228*(3*x(i)/L - 2.1)); % [kg/m3] - Density
        T(1,i) = T_0*(0.93968 + 0.0622*(3*x(i)/L - 2.1)); % [K] - Temperature

    end
end
u(1,:) = 0.59*rho_0*At*a_0./(rho(1,:).*A); % [m/s] - Velocity
e(1,:) = cv_0*T(1,:); % [m2/s2] - Internal specific energy
p(1,:) = rho(1,:).*Rg_0.*T(1,:); % [Pa] - Pressure

%% 4 - Solver

% Preparing flow-field variables structure:
vars.rho = rho; % [kg/m3] - Density
vars.u = u; % [m/s] - Velocity
vars.e = e; % [m2/s2] - Specific energy
vars.T = T; % [K] - Temperature
vars.p = p; % [Pa] - Pressure

% Solving:
tic
switch solver
    case "Godunov"
        [vars,gas,tvec] = SOLVER_Godunov(vars,A_f,dAdx_f,x,dx,t_end,sigma,rho_0,NaN,p_0,p_out,BCs,Riemann_solver, ...
            gas_model,Rg_0,gamma_0,species,n_0);
    case "MUSCL_Hancock"
        [vars,gas,tvec] = SOLVER_MUSCL_Hancock(vars,A_f,dAdx_f,x,dx,t_end,sigma,w,rho_0,NaN,p_0,p_out,BCs,Riemann_solver, ...
            gas_model,Rg_0,gamma_0,species,n_0);
end
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

%% 5 - Postprocess

% Compute other interesting data (from numerical solution):
a = sqrt(gamma.*Rg.*T); % [m/s] - Sound speed
M = u./a; % [] - Mach number
mdot = rho.*A.*u; % [kg/s] - Mass flow
T_0_num = T./M2Tr(M,gamma); % [K] - Stagnation temperature
p_0_num = p./M2pr(M,gamma); % [Pa] - Stagnation pressure

%% 6 - Exact analytical solution

% Prepare input data to generate the exact solution:
N_exact = 1000; % Number of data points
x_exact = linspace(0,L,N_exact); % [m] - Vector of grid coordinates
A_exact = arrayfun(A_f,x_exact); % [m2] - Cross-section area

% Call the analytical solution for perfect gases considering the presence
% of a shock wave inside the nozzle:
[M_exact,p_exact,rho_exact,T_exact,u_exact,p_0_exact,rho_0_exact,T_0_exact] = Nozzle_shock_solution(A_f,dAdx_f,At, ...
    x_exact,Rg_0,gamma_0,p_0,rho_0,T_0,p_out);

% Exact mass flow rate:
mdot_exact = rho_exact.*A_exact'.*u_exact; % [kg/s] - Mass flow

%% 7 - Read Anderson's numerical results

% Read text file with numerical results:
Anderson = readtable('RESULTS_Shock_nozzle_Anderson.txt','Delimiter','\t');

%% 8 - Plotting variable distributions

figure

% Temperature:
subplot(2,2,1)
hold on
plot(x,T(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,T_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
plot(Anderson.x_L*L/3,Anderson.T_T0*T_0,"-o","DisplayName","Anderson results","Color","red")
hold off
xlabel("x [m]")
ylabel("T [K]")
grid on
legend("Location","best")
title("Temperature")

% Pressure:
subplot(2,2,2)
hold on
plot(x,p(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,p_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
plot(Anderson.x_L*L/3,Anderson.p_p0*p_0,"-o","DisplayName","Anderson results","Color","red")
hold off
xlabel("x [m]")
ylabel("p [Pa]")
grid on
legend("Location","best")
title("Pressure")

% Velocity:
subplot(2,2,3)
hold on
plot(x,u(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,u_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
plot(Anderson.x_L*L/3,Anderson.V_a0*a_0,"-o","DisplayName","Anderson results","Color","red")
hold off
xlabel("x [m]")
ylabel("u [m/s]")
grid on
legend("Location","best")
title("Velocity")

% Mach number:
subplot(2,2,4)
hold on
plot(x,M(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,M_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
plot(Anderson.x_L*L/3,Anderson.M,"-o","DisplayName","Anderson results","Color","red")
hold off
xlabel("x [m]")
ylabel("Mach")
grid on
legend("Location","best")
title("Mach number")

%% 9 - Plotting mass flow and stagnation variables

figure

% Mass flow:
subplot(2,2,[1,2])
hold on
plot(x,mdot(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,mdot_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
plot(Anderson.x_L*L/3,Anderson.mdot*(rho_0*At*a_0),"-o","DisplayName","Anderson results","Color","red")
hold off
xlabel("x [m]")
ylabel("Mass flow [kg/s]")
grid on
legend("Location","best")
title("Mass flow")

% Stagnation temperature:
subplot(2,2,3)
hold on
plot(x,T_0_num(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,T_0_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("x [m]")
ylabel("T_0 [K]")
grid on
legend("Location","best")
title("Stagnation temperature")

% Stagnation pressure:
subplot(2,2,4)
hold on
plot(x,p_0_num(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,p_0_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("x [m]")
ylabel("p_0 [Pa]")
grid on
legend("Location","best")
title("Stagnation pressure")

% ----------------------------------------------------------------------- %