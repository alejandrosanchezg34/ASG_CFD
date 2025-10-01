% ----------------------------------------------------------------------- %
% Test case: Prandtl-Meyer expansion
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
solver = "MUSCL_Hancock"; % Either "Godunov" or "MUSCL_Hancock"

% Select Riemann solver:
Riemann_solver = "HLLC"; % Either "HLLC" or "exact"

% Configuration:
% - Supersonic inflow: fix all inlet variables (density, velocity and
% pressure).
% - Unconstrained outflow: let all variables float (density, velocity and
% pressure).
BCs = "supersonic_in_free_out";

% Numerical data:
N = 100; % Number of grid nodes
sigma = 0.9; % Stability criteria (CFL)
t_end = 10E-2; % [s] - Desired final time
w = 0; % [] - Constant parameter for MUSCL-Hancock's slope vectors (can go from -1 to 1)

% Geometric data:
fan_angle = 25; % [deg] - PM expansion fan angle
x_fan = 2; % [m] - x coordinate where the divergent region starts
A_i = 1; % [m2] - Inlet area
L = 10; % [m] - Domain length

% Gas properties:
gas_model = "constant"; % Constant gas model: properties and composition do not depend on T
species = NaN; % Gas species (NOT NEEDED for the constant gas model)
n_inf = NaN; % Gas composition (NOT NEEDED for the constant gas model)
Rg_inf = 287; % [J/(kgK)] - Gas constant
gamma_inf = 1.4; % [] - Gas adiabatic coefficient
cv_inf = getcv(Rg_inf,gamma_inf); % [J/(kgK)] - Gas specific heat at constant volume

% Flow properties:
rho_inf = 1.225; % [kg/m3] - Free stream density
p_inf = 1E5; % [Pa] - Free stream pressure
M_inf = 2; % [] - Free stream Mach number
T_inf = p_inf/(rho_inf*Rg_inf); % [K] - Free stream temperature
e_inf = cv_inf*T_inf; % [m2/s2] - Free stream specific energy
a_inf = sqrt(gamma_inf*Rg_inf*T_inf); % [m/s] - Free stream sound speed
u_inf = M_inf*a_inf; % [m/s] - Free stream velocity

%% 2 - Exact solution

% Prandtl-Meyer expansion given turn angle:
N_exact = 1000; % Number of data points
[M_exact,p_exact,rho_exact,T_exact,u_exact,ang_exact] = Prandtl_Meyer_solution(p_inf,T_inf,M_inf,Rg_inf,gamma_inf,fan_angle,N_exact);

%% 3 - 1D mesh generation

% Generating mesh:
[x,dx] = generate_mesh(L,N,"Godunov"); % [m], [m]

% Required inlet and exit areas:
Ai2At = M2Ar(M_inf,gamma_inf); % [] - Inlet to throat area ratio
Ae2At = M2Ar(M_exact(end),gamma_inf); % [] - Exit to throat area ratio
Ae2Ai = Ae2At/Ai2At; % [] - Exit to inlet area ratio
A_e = A_i*Ae2Ai; % [m2] - Exit area

% Generating cross-section area distribution:
A_f = @(x) (x <= x_fan)*A_i + ...
    (x > x_fan)*((A_e - A_i)/(L - x_fan)*(x - x_fan) + A_i); % [m2] - Cross-section area analytical function
dAdx_f = @(x) (x <= x_fan)*(0) + ...
    (x > x_fan)*((A_e - A_i)/(L - x_fan)); % Cross-section area derivative analytical function

A = arrayfun(A_f,x); % [m2] - Cross-section area distribution vector

%% 4 - Initial conditions

% Initializing variable matrices:
rho = zeros(1,N);
u = zeros(1,N);
p = zeros(1,N);
T = zeros(1,N);
e = zeros(1,N);

% Initial conditions (at t = 0):
rho(1,:) = rho_inf; % [kg/m3] - Density
u(1,:) = u_inf; % [m/s] - Velocity
p(1,:) = p_inf; % [Pa] - Pressure
T(1,:) = T_inf; % [K] - Temperature
e(1,:) = e_inf; % [m2/s2] - Internal specific energy

%% 5 - Solver

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
        [vars,gas,tvec] = SOLVER_Godunov(vars,A_f,dAdx_f,x,dx,t_end,sigma,rho_inf,u_inf,p_inf,NaN,BCs,Riemann_solver, ...
            gas_model,Rg_inf,gamma_inf,species,n_inf);
    case "MUSCL_Hancock"
        [vars,gas,tvec] = SOLVER_MUSCL_Hancock(vars,A_f,dAdx_f,x,dx,t_end,sigma,w,rho_inf,u_inf,p_inf,NaN,BCs,Riemann_solver, ...
            gas_model,Rg_inf,gamma_inf,species,n_inf);
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

%% 6 - Postprocess

% Compute other data from numerical solution:
a = sqrt(gamma.*Rg.*T); % [m/s] - Sound speed
M = u./a; % [] - Mach number
theta = rad2deg(PM_func(M,gamma) - PM_func(M_inf,gamma_inf)); % [deg] - Equivalent deflection angle

%% 7 - Plotting the resulting numerical distributions

figure

% Temperature at last time step:
subplot(2,2,1)
hold on
plot(x,T(end,:),"-o","DisplayName","Numerical results")
yline(T_exact(end),"--","DisplayName","Exact solution")
hold off
xlabel("x [m]")
ylabel("T [K]")
grid on
legend("Location","best")
title("Temperature")

% Pressure at last time step:
subplot(2,2,2)
hold on
plot(x,p(end,:),"-o","DisplayName","Numerical results")
yline(p_exact(end),"--","DisplayName","Exact solution")
hold off
xlabel("x [m]")
ylabel("p [Pa]")
grid on
legend("Location","best")
title("Pressure")

% Mach number at last time step:
subplot(2,2,3)
hold on
plot(x,M(end,:),"-o","DisplayName","Numerical results")
yline(M_exact(end),"--","DisplayName","Exact solution")
hold off
xlabel("x [m]")
ylabel("Mach")
grid on
legend("Location","best")
title("Mach number")

% Equivalent deflection angle at last time step:
subplot(2,2,4)
hold on
plot(x,theta(end,:),"-o","DisplayName","Numerical results")
yline(ang_exact(end),"--","DisplayName","Exact solution")
hold off
xlabel("x [m]")
ylabel("\theta [deg]")
grid on
legend("Location","best")
title("Equivalent turning angle")

%% 8 - Plotting the variables as functions of Mach

figure

% Pressure over Mach:
subplot(2,2,[1 2])
hold on
plot(M(end,:),p(end,:),"-o","DisplayName","Numerical results")
plot(M_exact,p_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("Mach")
ylabel("p [Pa]")
grid on
legend("Location","best")
title("Pressure")

% Temperature over Mach:
subplot(2,2,3)
hold on
plot(M(end,:),T(end,:),"-o","DisplayName","Numerical results")
plot(M_exact,T_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("Mach")
ylabel("T [K]")
grid on
legend("Location","best")
title("Temperature")

% Density over Mach:
subplot(2,2,4)
hold on
plot(M(end,:),rho(end,:),"-o","DisplayName","Numerical results")
plot(M_exact,rho_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("Mach")
ylabel("\rho [kg/m^3]")
grid on
legend("Location","best")
title("Density")

% ----------------------------------------------------------------------- %