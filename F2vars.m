% ----------------------------------------------------------------------- %
% Function that, given the flux vector F, returns the corresponding
% flow-field variables. Note that a non-linear system of equations is
% solved.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - F: flux vector ([3 x 1] or [3 x N]).
% - A: cross-section area [m2] (may be given as a [1 x N] vector).
% - gas_model: string variable that specifies the gas model selected. It
% can be either "constant", "frozen_flow" or "shifting_flow".
% - Rg_ct: gas constant to be used by the constant gas model [J/(kgK)].
% - gamma_ct: gas adiabatic coefficient to be used by the constant gas
% model [].
% - species: cell array with the names (strings) of the species involved in
% the gas. It is not needed for the "constant" gas model.
% - n_0: initial/chamber gas composition [mols] (given as a [species x 1]
% vector).
% - vars_ini: initial guesses for the variables, consisting of a vector of
% the form [rho; u; p; e; T] ([5 x 1] or [5 x N]).
% ----------------------------------------------------------------------- %
% Output:
% - rho: density [kg/m3].
% - u: velocity [m/s].
% - e: internal specific energy [m2/s2].
% - T: temperature [K].
% - p: pressure [Pa].
% * May be given as vectors [1 x N].
% ----------------------------------------------------------------------- %

function [rho,u,e,T,p] = F2vars(F,A,gas_model,Rg_ct,gamma_ct,species,n_0,vars_ini)

% Number of grid points:
N = size(F,2);

% Initializing variable vectors:
rho = zeros(1,N);
u = zeros(1,N);
e = zeros(1,N);
T = zeros(1,N);
p = zeros(1,N);

% Solve component-by-component:
for i = 1:N

    % Variables of the system of equations to solve:
    % z(1) = rho;
    % z(2) = u;
    % z(3) = p;
    % z(4) = e;
    % z(5) = T;

    % Define gas properties as functions of both z(3) = p and z(5) = T:
    Rg = @(z) getRg(gas_model,Rg_ct,species,n_0,z(5),z(3)); % [J/(kgK)] - Gas constant
    gamma = @(z) getgamma(gas_model,gamma_ct,species,n_0,z(5),z(3)); % [] - Gas adiabatic coefficient
    cv = @(z) getcv(Rg(z),gamma(z)); % [J/(kgK)] - Gas specific heat at constant volume

    % Define system of equations to solve:
    Func = @(z) [z(1)*A(1,i)*z(2) - F(1,i);
                 z(1)*A(1,i)*z(2)^2 + z(3)*A(1,i) - F(2,i);
                 z(1)*(z(4) + z(2)^2/2)*z(2)*A(1,i) + z(3)*A(1,i)*z(2) - F(3,i);
                 z(4) - cv(z)*z(5);
                 z(3) - state_equation(z(1),Rg(z),z(5))];

    % Solve system:
    z = fsolve(Func,vars_ini(:,i),optimoptions('fsolve','Display','off'));
    rho(1,i) = z(1); % [kg/m3] - Density
    u(1,i) = z(2); % [m/s] - Velocity
    p(1,i) = z(3); % [Pa] - Pressure
    e(1,i) = z(4); % [m2/s2] - Internal specific energy
    T(1,i) = z(5); % [K] - Temperature

end

end

% ----------------------------------------------------------------------- %