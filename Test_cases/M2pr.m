% ----------------------------------------------------------------------- %
% This function provides the pressure ratio p/p0 associated to a given Mach
% number M, where p0 is the stagnation pressure.
% ----------------------------------------------------------------------- %
% Input:
% - M: Mach number [].
% - gamma: gas adiabatic coefficient [].
% ----------------------------------------------------------------------- %
% Output:
% - p_p0: corresponding p/p0 ratio [].
% ----------------------------------------------------------------------- %

function p_p0 = M2pr(M,gamma)

p_p0 = (1 + (gamma-1)./2.*M.^2).^(-gamma./(gamma-1));

end

% ----------------------------------------------------------------------- %