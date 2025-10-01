% ----------------------------------------------------------------------- %
% This function provides the result of the Prandtl-Meyer function for a
% given Mach number and gas adiabatic coefficient.
% ----------------------------------------------------------------------- %
% Input:
% - M: Mach number [].
% - gamma: gas adiabatic coefficient [].
% ----------------------------------------------------------------------- %
% Output:
% - nu: Prandtl-Meyer function result for the given Mach number M [rad].
% ----------------------------------------------------------------------- %

function nu = PM_func(M,gamma)

nu = sqrt((gamma+1)./(gamma-1)).*atan(sqrt((gamma-1)./(gamma+1).*(M.^2 - 1))) - atan(sqrt(M.^2 - 1));

end

% ----------------------------------------------------------------------- %