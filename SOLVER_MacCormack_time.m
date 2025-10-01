% ----------------------------------------------------------------------- %
% This function implements the time-marching MacCormack method in order to
% obtain the time-marching solution of the quasi-1D Euler equations. If a
% constant area distribution is given, the pure-1D Euler equations will be
% solved, instead.
% ----------------------------------------------------------------------- %
% Input:
% - vars: structure containing the arrays with the flow-field variables.
% - A_f: cross-section area analytical function (function handle).
% - dAdx_f: cross-section area derivative analytical function (function handle).
% - x: grid nodes coordinates vector [1 x N].
% - dx: grid's cell width [m].
% - t_end: final simulation time [s].
% - sigma: stability parameter.
% - rho_in: if needed, inlet density to be fixed [kg/m3].
% - u_in: if needed, inlet velocity to be fixed [m/s].
% - p_in: if needed, inlet pressure to be fixed [Pa].
% - p_out: if needed, outlet pressure to be fixed [Pa].
% - Cx: artificial viscosity coefficient.
% - BCs: string specifying the type of boundary conditions to impose.
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

function [vars,gas,tvec] = SOLVER_MacCormack_time(vars,A_f,dAdx_f,x,dx,t_end,sigma, ...
    rho_in,u_in,p_in,p_out,Cx,BCs,gas_model,Rg_ct,gamma_ct,species,n_0)

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

% Build 'A' and 'dAdx' vectors:
A = arrayfun(A_f,x); % [m2] - Cross-section area
dAdx = arrayfun(dAdx_f,x); % Cross-section area derivatives vector

% Initializing solver data and arrays:
j = 1; % Initializing time steps
tvec = 0; % Initializing time vector
U = zeros(3,N,1); % Initializing solution U matrix
F = zeros(3,N,1); % Initializing flux F matrix
J = zeros(3,N,1); % Initializing source terms J matrix

% Transient (time-marching) solver algorithm:
while (tvec(end) < t_end && j < Nit)

    % Initializing vectors and matrices:
    dU_dt_now = zeros(3,N); % Current time derivatives
    dU_dt_pred = zeros(3,N); % Predicted time derivatives
    U_pred = zeros(3,N); % Predicted solution vectors
    S_now = zeros(3,N); % Current artificial viscosity
    S_pred = zeros(3,N); % Predicted artificial viscosity

    % 1- Compute U from vars:
    U(:,:,j) = vars2U(rho(j,:),u(j,:),e(j,:),A);

    % 2- Compute F from vars:
    F(:,:,j) = vars2F(rho(j,:),u(j,:),e(j,:),p(j,:),A);

    % 3- Compute J from vars:
    J(:,:,j) = vars2J(p(j,:),dAdx);

    % 4- Compute time step required for stability:
    dt = stability_dt_Sod(u(j,:),T(j,:),Rg(j,:),gamma(j,:),sigma,dx); % [s]
    
    % Correct 'dt' for the last time step:
    if (tvec(j) + dt > t_end)
        dt = t_end - tvec(j); % [s]
    end

    % 5- Employ MacCormack time-marching scheme:
    
    % --- Predictor step
    
    % a) Time derivatives at current instant:
    for i = 2:(N-1)
        dU_dt_now(:,i) = -(1/dx)*(F(:,i+1,j) - F(:,i,j)) + J(:,i,j);
    end
    
    % b) Artificial viscosity for the predictor step:
    for i = 2:(N-1)
        S_now(:,i) = art_visc_MC_1D(p(j,:),U(:,:,j),Cx,i);
    end
    
    % c) Predicted value of the solution vector U:
    for i = 2:(N-1)
        U_pred(:,i) = U(:,i,j) + dU_dt_now(:,i)*dt + S_now(:,i);
    end

    % Obtain predicted flow-field variables (only internal nodes):
    [rho_pred,u_pred,e_pred,T_pred,p_pred] = U2vars(U_pred,A,Rg(j,:),gamma(j,:));

    % Enforce boundary conditions on predicted variables:
    [rho_pred,u_pred,e_pred,~,p_pred] = BCs_MacCormack(rho_pred,u_pred, ...
        e_pred,T_pred,p_pred,Rg(j,:),gamma(j,:),rho_in,u_in,p_in,p_out,BCs);

    % Compute U_pred from predicted variables:
    U_pred = vars2U(rho_pred,u_pred,e_pred,A);

    % Compute F_pred from predicted variables:
    F_pred = vars2F(rho_pred,u_pred,e_pred,p_pred,A);

    % Compute J_pred from predicted variables:
    J_pred = vars2J(p_pred,dAdx);

    % --- Corrector step
    
    % a) Predicted time derivative:
    for i = 2:(N-1)
        dU_dt_pred(:,i) = -(1/dx)*(F_pred(:,i) - F_pred(:,i-1)) + J_pred(:,i);
    end
    
    % b) Artificial viscosity for the corrector step:
    for i = 2:(N-1)
        S_pred(:,i) = art_visc_MC_1D(p_pred,U_pred,Cx,i);
    end
    
    % c) Average time derivative:
    dU_dt_avg = (1/2)*(dU_dt_now + dU_dt_pred);
    
    % --- Final solution
    for i = 2:(N-1)
        U(:,i,j+1) = U(:,i,j) + dU_dt_avg(:,i)*dt + S_pred(:,i);
    end

    % Obtain corrected (new step) flow-field variables (only internal nodes):
    [rho(j+1,:),u(j+1,:),e(j+1,:),T(j+1,:),p(j+1,:)] = U2vars(U(:,:,j+1),A,Rg(j,:),gamma(j,:));

    % Enforce boundary conditions on corrected (new step) variables:
    [rho(j+1,:),u(j+1,:),e(j+1,:),T(j+1,:),p(j+1,:)] = BCs_MacCormack(rho(j+1,:),u(j+1,:), ...
        e(j+1,:),T(j+1,:),p(j+1,:),Rg(j,:),gamma(j,:),rho_in,u_in,p_in,p_out,BCs);

    % 6- Obtain gas properties:
    Rg(j+1,:) = getRg(gas_model,Rg_ct,species,n_0,T(j+1,:),p(j+1,:)); % [J/(kgK)] - Gas constant
    gamma(j+1,:) = getgamma(gas_model,gamma_ct,species,n_0,T(j+1,:),p(j+1,:)); % [] - Gas adiabatic coefficient

    % 7- Update time vector:
    tvec = [tvec, tvec(j) + dt];
    fprintf("Current time step: %.4e [s] | Simulated time: %.4e [s] \n",dt,tvec(end));
    
    % 8- Go for next time instant:
    j = j + 1;

    % 9- Check convergence to the steady state:
    if check_convergence(rho)
        disp("Steady state has been reached. Stopping calculations...");
        break;
    end

end

% Saving flow-field variables:
vars.rho = rho; % [kg/m3] - Density
vars.u = u; % [m/s] - Velocity
vars.e = e; % [m2/s2] - Specific energy
vars.T = T; % [K] - Temperature
vars.p = p; % [Pa] - Pressure

% Saving gas properties:
gas.Rg = Rg; % [J/(kgK)] - Gas constant
gas.gamma = gamma; % [] - Gas adiabatic coefficient

end

% ----------------------------------------------------------------------- %