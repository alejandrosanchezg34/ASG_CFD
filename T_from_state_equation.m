% ----------------------------------------------------------------------- %
% Function that applies the state equation for ideal gases in order to get
% the temperature, for a given density and pressure.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density [kg/m3] (may be given as a [1 x N] vector).
% - Rg: gas constant [J/(kgK)] (may be given as a [1 x N] vector).
% - p: pressure [Pa] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %
% Output:
% - T: temperature [K] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %

function T = T_from_state_equation(rho,Rg,p)

% State equation for ideal gases:
T = p./(rho.*Rg); % [K] - Temperature

end

% ----------------------------------------------------------------------- %