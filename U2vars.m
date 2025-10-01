% ----------------------------------------------------------------------- %
% Function that, given the solution vector U, returns the corresponding
% flow-field variables.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - U: solution vector ([3 x 1] or [3 x N]).
% - A: cross-section area [m2] (may be given as a [1 x N] vector).
% - Rg: gas constant [J/(kgK)] (may be given as a [1 x N] vector).
% - gamma: gas adiabatic coefficient [] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %
% Output:
% - rho: density [kg/m3].
% - u: velocity [m/s].
% - e: internal specific energy [m2/s2].
% - T: temperature [K].
% - p: pressure [Pa].
% * May be given as vectors [1 x N].
% ----------------------------------------------------------------------- %

function [rho,u,e,T,p] = U2vars(U,A,Rg,gamma)

% --- Compute flow-field variables available from U
rho = U(1,:)./A(1,:); % [kg/m3] - Density
u = U(2,:)./U(1,:); % [m/s] - Velocity
e = U(3,:)./U(1,:) - u.^2./2; % [m2/s2] - Internal specific energy

% --- Compute temperature T from specific energy e
cv = getcv(Rg,gamma); % [J/(kgK)] - Gas specific heat at constant volume
T = e./cv; % [K] - Temperature

% --- Compute pressure p from the equation of state
p = state_equation(rho,Rg,T); % [Pa] - Pressure

end

% ----------------------------------------------------------------------- %