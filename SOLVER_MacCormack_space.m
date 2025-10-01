% ----------------------------------------------------------------------- %
% This function implements the space-marching MacCormack method in order to
% obtain the steady solution of the quasi-1D Euler equations. If a constant
% area distribution is given, the pure-1D Euler equations will be solved,
% instead.
% This scheme can only be applied to locally supersonic flows everywhere.
% ----------------------------------------------------------------------- %
% Input:
% - vars: structure containing the arrays with the flow-field variables.
% - A_f: cross-section area analytical function (function handle).
% - dAdx_f: cross-section area derivative analytical function (function handle).
% - x: grid nodes coordinates vector [1 x N].
% - dx: grid's cell width [m].
% - gas_model: string variable that specifies the gas model selected. It
% can be either "constant", "frozen_flow" or "shifting_flow".
% - Rg_ct: gas constant to be used by the constant gas model [J/(kgK)].
% - gamma_ct: gas adiabatic coefficient to be used by the constant gas
% model [].
% - species: cell array with the names (strings) of the species involved in
% the gas. It is not needed for the "constant" gas model.
% - n_0: initial/chamber gas composition [mols] (given as a [species x 1]
% vector).
% ----------------------------------------------------------------------- %
% Output:
% - vars: structure containing the arrays with the flow-field variables.
% - gas: structure containing the arrays with the gas properties.
% ----------------------------------------------------------------------- %

function [vars,gas] = SOLVER_MacCormack_space(vars,A_f,dAdx_f,x,dx,gas_model,Rg_ct,gamma_ct,species,n_0)

% Decoding variables and imposing initial conditions:
rho = vars.rho; % [kg/m3] - Density
u = vars.u; % [m/s] - Velocity
e = vars.e; % [m2/s2] - Specific energy
T = vars.T; % [K] - Temperature
p = vars.p; % [Pa] - Pressure

% Number of nodes being evaluated:
N = size(rho,2);

% Build 'A' and 'dAdx' vectors:
A = arrayfun(A_f,x); % [m2] - Cross-section area
dAdx = arrayfun(dAdx_f,x); % Cross-section area derivatives vector

% Initializing solver data and arrays:
F = zeros(3,N); % Initializing flux F matrix
J = zeros(3,N); % Initializing source terms J matrix

% Impose inlet data:
F(:,1) = vars2F(rho(1,1),u(1,1),e(1,1),p(1,1),A(1,1));

% Steady (space-marching) solver algorithm:
for i = 1:(N-1)

    % 1- Compute J from vars:
    J(:,i) = vars2J(p(1,i),dAdx(1,i));

    % 2- Employ MacCormack space-marching scheme:
    
    % --- Predictor step
    
    % a) Spatial derivative at current node 'i':
    dF_dx_now = J(:,i);

    % b) Predicted value of the flux vector F:
    F_pred = F(:,i) + dF_dx_now*dx;

    % Obtain predicted flow-field variables:
    vars_old = [rho(1,i);u(1,i);p(1,i);e(1,i);T(1,i)]; % Initial guesses for the variables
    [rho_pred,u_pred,e_pred,T_pred,p_pred] = F2vars(F_pred,A(1,i+1),gas_model,Rg_ct,gamma_ct,species,n_0,vars_old);

    % Compute J_pred from predicted variables:
    J_pred = vars2J(p_pred,dAdx(1,i+1));

    % --- Corrector step
    
    % a) Predicted spatial derivative:
    dF_dx_pred = J_pred;

    % b) Average spatial derivative:
    dF_dx_avg = (1/2)*(dF_dx_now + dF_dx_pred);

    % --- Final solution
    F(:,i+1) = F(:,i) + dF_dx_avg*dx;

    % Obtain corrected (new node) flow-field variables:
    vars_old = [rho_pred;u_pred;p_pred;e_pred;T_pred]; % Initial guesses for the variables
    [rho(1,i+1),u(1,i+1),e(1,i+1),T(1,i+1),p(1,i+1)] = F2vars(F(:,i+1),A(1,i+1),gas_model,Rg_ct,gamma_ct,species,n_0,vars_old);

    fprintf("Solved node %d of %d \n",i+1,N);

end

% Saving flow-field variables:
vars.rho = rho; % [kg/m3] - Density
vars.u = u; % [m/s] - Velocity
vars.e = e; % [m2/s2] - Specific energy
vars.T = T; % [K] - Temperature
vars.p = p; % [Pa] - Pressure

% Saving gas properties:
gas.Rg = getRg(gas_model,Rg_ct,species,n_0,T,p); % [J/(kgK)] - Gas constant
gas.gamma = getgamma(gas_model,gamma_ct,species,n_0,T,p); % [] - Gas adiabatic coefficient

end

% ----------------------------------------------------------------------- %