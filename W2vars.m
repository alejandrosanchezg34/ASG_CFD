% ----------------------------------------------------------------------- %
% Function that, given the Riemann problem solution vector W, returns the
% corresponding flow-field variables.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - W: solution vector ([3 x 1] or [3 x N]).
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

function [rho,u,e,T,p] = W2vars(W,Rg,gamma)

% Compute flow-field variables from 'W':
rho = W(1,:); % [kg/m3] - Density
u = W(2,:); % [m/s] - Velocity
p = W(3,:); % [Pa] - Pressure

% Get gas specific heat at constant volume:
cv = getcv(Rg,gamma); % [J/(kgK)]

% Compute corresponding temperature T:
T = T_from_state_equation(rho,Rg,p); % [K] - Temperature

% Compute corresponding internal specific energy e:
e = cv.*T; % [m2/s2] - Internal specific energy

end

% ----------------------------------------------------------------------- %