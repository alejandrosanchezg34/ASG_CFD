% ----------------------------------------------------------------------- %
% This function provides the density ratio rho/rho0 associated to a given
% Mach number M, where rho0 is the stagnation density.
% ----------------------------------------------------------------------- %
% Input:
% - M: Mach number [].
% - gamma: gas adiabatic coefficient [].
% ----------------------------------------------------------------------- %
% Output:
% - rho_rho0: corresponding rho/rho0 ratio [].
% ----------------------------------------------------------------------- %

function rho_rho0 = M2rhor(M,gamma)

rho_rho0 = (1 + (gamma-1)./2.*M.^2).^(-1./(gamma-1));

end

% ----------------------------------------------------------------------- %