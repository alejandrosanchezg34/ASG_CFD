% ----------------------------------------------------------------------- %
% Test case: Sod shock tube problem
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

% Configuration:
% - Transmissive inlet: impose zero gradient.
% - Transmissive outlet: impose zero gradient.
BCs = "zero_gradient";

% Numerical data:
Cx = 0.3; % [] - Artificial viscosity coefficient (for MacCormack's scheme)
N = 101; % Number of grid nodes
sigma = 0.9; % Stability criteria (CFL)
t_end = 200E-3; % [s] - Desired final time

% Geometric data:
L = 1; % [m] - Tube length

% Gas properties:
gas_model = "constant"; % Constant gas model: properties and composition do not depend on T
species = NaN; % Gas species (NOT NEEDED for the constant gas model)
n_0 = NaN; % Gas composition (NOT NEEDED for the constant gas model)
Rg_0 = 287; % [J/(kgK)] - Gas constant
gamma_0 = 1.4; % [] - Gas adiabatic coefficient
cv_0 = getcv(Rg_0,gamma_0); % [J/(kgK)] - Gas specific heat at constant volume

% Flow properties:
rhoL = 1.0; % [kg/m3] - Left region density
pL = 1.0; % [Pa] - Left region pressure
uL = 0.0; % [m/s] - Left region velocity
WL = [rhoL; uL; pL]; % Left region data

rhoR = 0.125; % [kg/m3] - Right region density
pR = 0.1; % [Pa] - Right region pressure
uR = 0.0; % [m/s] - Right region velocity
WR = [rhoR; uR; pR]; % Right region data

%% 2 - 1D mesh generation

% Generating mesh:
[x,dx] = generate_mesh(L,N,"MacCormack",-L/2); % [m], [m]

% Generating tube cross-section area:
A_f = @(x) 1; % [m2] - Cross-section area analytical function
dAdx_f = @(x) 0; % Cross-section area derivative analytical function

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
    if x(i) <= 0 % Left side of the membrane
        rho(1,i) = rhoL; % [kg/m3]
        u(1,i) = uL; % [m/s]
        p(1,i) = pL; % [Pa]
    else % Right side of the membrane
        rho(1,i) = rhoR; % [kg/m3]
        u(1,i) = uR; % [m/s]
        p(1,i) = pR; % [Pa]
    end
end
T(1,:) = p(1,:)./(rho(1,:)*Rg_0); % [K]
e(1,:) = cv_0*T(1,:); % [m2/s2]

%% 4 - Solver

% Preparing flow-field variables structure:
vars.rho = rho; % [kg/m3] - Density
vars.u = u; % [m/s] - Velocity
vars.e = e; % [m2/s2] - Specific energy
vars.T = T; % [K] - Temperature
vars.p = p; % [Pa] - Pressure

% Solving:
tic
[vars,gas,tvec] = SOLVER_MacCormack_time(vars,A_f,dAdx_f,x,dx,t_end,sigma,NaN,NaN,NaN,NaN,Cx,BCs, ...
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

%% 5 - Exact solution

% Prepare input data to generate the exact solution:
N_exact = 1000; % Number of data points
x_exact = linspace(-L/2,L/2,N_exact); % [m] - Vector of grid coordinates

% Calling exact solution to the Sod shock tube problem:
WSol = zeros(3,N_exact); % Initializing solution vectors
for i = 1:N_exact
    WSol(:,i) = Riemann_solution(x_exact(i),tvec(end),gamma_0,WL,WR);
end

% Decode exact flow-field variables:
[rho_exact,u_exact,e_exact,T_exact,p_exact] = W2vars(WSol,Rg_0,gamma_0);

%% 6 - Plotting results

figure

% Density at last time step:
subplot(2,2,1)
hold on
plot(x,rho(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,rho_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("x [m]")
ylabel("\rho [kg/m^3]")
grid on
legend("Location","best")
title("Density")
xticks([-0.5 -0.25 0 0.25 0.5])

% Velocity at last time step:
subplot(2,2,2)
hold on
plot(x,u(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,u_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("x [m]")
ylabel("u [m/s]")
grid on
legend("Location","best")
title("Velocity")
xticks([-0.5 -0.25 0 0.25 0.5])

% Pressure at last time step:
subplot(2,2,3)
hold on
plot(x,p(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,p_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("x [m]")
ylabel("p [Pa]")
grid on
legend("Location","best")
title("Pressure")
xticks([-0.5 -0.25 0 0.25 0.5])

% Internal specific energy at last time step:
subplot(2,2,4)
hold on
plot(x,e(end,:),"-o","DisplayName","Numerical results")
plot(x_exact,e_exact,"LineWidth",1,"DisplayName","Exact solution","Color","black")
hold off
xlabel("x [m]")
ylabel("e [m^2/s^2]")
grid on
legend("Location","best")
title("Internal specific energy")
xticks([-0.5 -0.25 0 0.25 0.5])

% ----------------------------------------------------------------------- %