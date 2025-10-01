% ----------------------------------------------------------------------- %
% This function implements the second-order MUSCL-Hancock scheme applied to
% the second version of Godunov's first-order upwind method to obtain the
% time-marching solution of the quasi-1D Euler equations. If a constant
% area distribution is given, the pure-1D Euler equations will be solved.
% The Riemann solver to be used can either be the exact one or the HLLC
% approximate one.
% ----------------------------------------------------------------------- %
% Input:
% - vars: structure containing the arrays with the flow-field variables.
% - A_f: cross-section area analytical function (function handle).
% - dAdx_f: cross-section area derivative analytical function (function handle).
% - x: grid nodes coordinates vector [1 x N].
% - dx: grid's cell width [m].
% - t_end: final simulation time [s].
% - sigma: stability parameter.
% - w: constant parameter to compute slope vectors (can go from -1 to 1).
% - rho_in: if needed, inlet density to be fixed [kg/m3].
% - u_in: if needed, inlet velocity to be fixed [m/s].
% - p_in: if needed, inlet pressure to be fixed [Pa].
% - p_out: if needed, outlet pressure to be fixed [Pa].
% - BCs: string specifying the type of boundary conditions to impose.
% - Riemann_solver: Riemann problem solver to be used ("HLLC" or "exact").
% - gas_model: string variable with the name of the gas model selected.
% - Rg_ct: gas constant to be used with the constant gas model [J/(kgK)].
% - gamma_ct: gas adiabatic coefficient to be used with the constant gas
% model [].
% - species: cell array with the names (strings) of the species involved in
% the gas.
% - n_0: initial/chamber gas composition [mols] (given as a [species x 1]
% vector).
% ----------------------------------------------------------------------- %
% Output:
% - vars: structure containing the arrays with the flow-field variables.
% - gas: structure containing the arrays with the gas properties.
% - tvec: vector containing all simulated times.
% ----------------------------------------------------------------------- %

function [vars,gas,tvec] = SOLVER_MUSCL_Hancock(vars,A_f,dAdx_f,x,dx,t_end,sigma,w, ...
    rho_in,u_in,p_in,p_out,BCs,Riemann_solver,gas_model,Rg_ct,gamma_ct,species,n_0)

% ----------------------------------------------------------------------- %
% Set number of maximum time steps:
Nit = 1E6;
% ----------------------------------------------------------------------- %

% Decoding variables and imposing initial conditions:
rho = vars.rho; % [kg/m3] - Density
u = vars.u; % [m/s] - Velocity
e = vars.e; % [m2/s2] - Specific energy
T = vars.T; % [K] - Temperature
p = vars.p; % [Pa] - Pressure

% Number of nodes being evaluated:
N = size(rho,2);

% Get initial gas properties:
Rg = getRg(gas_model,Rg_ct,species,n_0,T,p); % [J/(kgK)] - Gas constant
gamma = getgamma(gas_model,gamma_ct,species,n_0,T,p); % [] - Gas adiabatic coefficient

% Adding additional nodes to evaluate BCs (nodes 1, 2 and N+3, N+4):
rho = [0, 0, rho, 0, 0]; % [kg/m3] - Density
u = [0, 0, u, 0, 0]; % [m/s] - Velocity
e = [0, 0, e, 0, 0]; % [m2/s2] - Specific energy
T = [0, 0, T, 0, 0]; % [K] - Temperature
p = [0, 0, p, 0, 0]; % [Pa] - Pressure

Rg = [Rg(1), Rg(1), Rg, Rg(end), Rg(end)]; % [J/(kgK)] - Gas constant
gamma = [gamma(1), gamma(1), gamma, gamma(end), gamma(end)]; % [] - Gas adiabatic coefficient

x = [x(1)-2*dx, x(1)-dx, x, x(end)+dx, x(end)+2*dx]; % [m] - Vector of grid coordinates

% Build 'A' and 'dAdx' vectors:
A = arrayfun(A_f,x); % [m2] - Cross-section area
dAdx = arrayfun(dAdx_f,x); % Cross-section area derivatives vector

% Initializing solver data and arrays:
j = 1; % Initializing time steps
tvec = 0; % Initializing time vector
U = zeros(3,N+4,1); % Initializing solution U matrix
F = zeros(3,N+4,1); % Initializing flux F matrix
J = zeros(3,N+4,1); % Initializing source terms J matrix

% Transient (time-marching) solver algorithm:
while (tvec(end) < t_end && j < Nit)

    % Approximate gas properties to those of the last time step (except for
    % the first time step):
    if j > 1
        Rg(j,:) = Rg(j-1,:); % [J/(kgK)] - Gas constant
        gamma(j,:) = gamma(j-1,:); % [] - Gas adiabatic coefficient
    end

    % 1- Apply boundary conditions:
    [rho(j,:),u(j,:),e(j,:),T(j,:),p(j,:)] = BCs_MUSCL(rho(j,:),u(j,:),e(j,:), ...
        T(j,:),p(j,:),Rg(j,:),gamma(j,:),rho_in,u_in,p_in,p_out,BCs);

    % 2- Obtain new gas properties:
    Rg(j,:) = getRg(gas_model,Rg_ct,species,n_0,T(j,:),p(j,:)); % [J/(kgK)] - Gas constant
    gamma(j,:) = getgamma(gas_model,gamma_ct,species,n_0,T(j,:),p(j,:)); % [] - Gas adiabatic coefficient

    % 3- Compute U from vars:
    U(:,:,j) = vars2U(rho(j,:),u(j,:),e(j,:),A);

    % 4- Compute F from vars:
    F(:,:,j) = vars2F(rho(j,:),u(j,:),e(j,:),p(j,:),A);

    % 5- Compute J from vars:
    J(:,:,j) = vars2J(p(j,:),dAdx);

    % 6- Find non-linear wave speeds throughout domain and Riemann problem
    % data:
    [SL,SR] = domain_speeds_Godunov(U(:,:,j),A,Rg(j,:),gamma(j,:),Riemann_solver);

    % 7- Time step required for stability:
    dt = stability_dt_waves(SL,SR,dx,sigma); % [s]
    
    % For last time step:
    if (tvec(j) + dt > t_end)
        dt = t_end - tvec(j); % [s]
    end

    % 8- Employ the time-marching MUSCL-Hancock scheme (based on Godunov's
    % upwind method):
    UiL_new = zeros(3,N+4); % Initializing
    UiR_new = zeros(3,N+4); % Initializing

    for i = 2:(N+3)

        % --- Reconstructing data and obtaining the boundary extrapolated
        % values for the internal nodes (and internal BC nodes):

        % Compute slope vector at cell 'i':
        di = slope_vector(i,U(:,:,j),w);

        % Compute boundary extrapolated values at cell 'i':
        UiL = U(:,i,j) - 0.5*di;
        UiR = U(:,i,j) + 0.5*di;

        % --- Evolve boundary extrapolated values by a time 1/2*dt:

        % Compute intercell cross-section areas:
        A_int_L = arrayfun(A_f,x(i) - 0.5*dx); % [m2] - Left intercell
        A_int_R = arrayfun(A_f,x(i) + 0.5*dx); % [m2] - Right intercell

        % Compute intercell gas constants:
        Rg_int_L = harmmean([Rg(j,i-1),Rg(j,i)]); % [J/(kgK)]
        Rg_int_R = harmmean([Rg(j,i),Rg(j,i+1)]); % [J/(kgK)]

        % Compute intercell gas adiabatic coefficients:
        gamma_int_L = harmmean([gamma(j,i-1),gamma(j,i)]); % []
        gamma_int_R = harmmean([gamma(j,i),gamma(j,i+1)]); % []

        % Fluxes associated to the boundary extrapolated values:
        F_UiL = U2F(UiL,A_int_L,Rg_int_L,gamma_int_L);
        F_UiR = U2F(UiR,A_int_R,Rg_int_R,gamma_int_R);

        % Evolve by a time 1/2*dt:
        UiL_new(:,i) = UiL + 0.5*dt/dx*(F_UiL - F_UiR) + 0.5*dt*J(:,i,j);
        UiR_new(:,i) = UiR + 0.5*dt/dx*(F_UiL - F_UiR) + 0.5*dt*J(:,i,j);
    end

    % 9- Solve Riemann problems to get intercell fluxes 'i+1/2':
    F_int = zeros(3,N+2); % Initializing

    for i = 2:(N+2)

        % Compute 'i+1/2' intercell variables:
        A_int = arrayfun(A_f,x(i) + 0.5*dx); % [m2] - Cross-section area
        Rg_int = harmmean([Rg(j,i),Rg(j,i+1)]); % [J/(kgK)] - Gas constant
        gamma_int = harmmean([gamma(j,i),gamma(j,i+1)]); % [] - Gas adiabatic coefficient

        % Get Riemann problem data:
        WL = U2W(UiR_new(:,i),A_int,Rg_int,gamma_int); % Left
        WR = U2W(UiL_new(:,i+1),A_int,Rg_int,gamma_int); % Right
        
        % Compute 'i+1/2' intercell flux:
        switch Riemann_solver
            case "HLLC"
                F_HLLC = HLLC_solver(WL,WR,Rg_int,gamma_int); % HLLC intercell flux
                F_int(:,i) = F_HLLC*A_int; % Convert to quasi-1D flux vector

            case "exact"
                W_int = Riemann_solution(0,dt,gamma_int,WL,WR); % Riemann problem solution
                F_int(:,i) = W2F(W_int,A_int,Rg_int,gamma_int); % Get flux vector
        end
    end

    % 10- Employ the second version of Godunov's method (time-marching
    % scheme) for internal nodes:
    for i = 3:(N+2)

        % Left and right intercell fluxes:
        Fint_L = F_int(:,i-1);
        Fint_R = F_int(:,i);

        % Next time instant:
        U(:,i,j+1) = U(:,i,j) + (dt/dx)*(Fint_L - Fint_R) + dt*J(:,i,j);

    end

    % 11- New step flow-field variables:
    [rho(j+1,:),u(j+1,:),e(j+1,:),T(j+1,:),p(j+1,:)] = U2vars(U(:,:,j+1),A,Rg(j,:),gamma(j,:));

    % 12- Update time vector:
    tvec = [tvec, tvec(j) + dt];
    fprintf("Current time step: %.4e [s] | Simulated time: %.4e [s] \n",dt,tvec(end));
    
    % 13- Go for next time step:
    j = j + 1;

    % 14- Check convergence to the steady state:
    if check_convergence(rho(:,3:(N+2)))
        disp("Steady state has been reached. Stopping calculations...");
        break;
    end

end

% Get gas properties at the last time step:
Rg(j,:) = getRg(gas_model,Rg_ct,species,n_0,T(j,:),p(j,:)); % [J/(kgK)] - Gas constant
gamma(j,:) = getgamma(gas_model,gamma_ct,species,n_0,T(j,:),p(j,:)); % [] - Gas adiabatic coefficient

% Saving flow-field variables (only internal nodes):
vars.rho = rho(:,3:(N+2)); % [kg/m3] - Density
vars.u = u(:,3:(N+2)); % [m/s] - Velocity
vars.e = e(:,3:(N+2)); % [m2/s2] - Specific energy
vars.T = T(:,3:(N+2)); % [K] - Temperature
vars.p = p(:,3:(N+2)); % [Pa] - Pressure

% Saving gas properties (only internal nodes):
gas.Rg = Rg(:,3:(N+2)); % [J/(kgK)] - Gas constant
gas.gamma = gamma(:,3:(N+2)); % [] - Gas adiabatic coefficient

end

% ----------------------------------------------------------------------- %