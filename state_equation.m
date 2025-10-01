% ----------------------------------------------------------------------- %
% Function that applies the state equation for ideal gases in order to get
% the pressure, for a given density and temperature.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density [kg/m3] (may be given as a [1 x N] vector).
% - Rg: gas constant [J/(kgK)] (may be given as a [1 x N] vector).
% - T: temperature [K] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %
% Output:
% - p: pressure [Pa] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %

function p = state_equation(rho,Rg,T)

% State equation for ideal gases:
p = rho.*Rg.*T; % [Pa] - Pressure

end

% ----------------------------------------------------------------------- %