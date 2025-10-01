% ----------------------------------------------------------------------- %
% This function provides the area ratio A/At associated to a given Mach
% number M, where At is the sonic throat area.
% ----------------------------------------------------------------------- %
% Input:
% - M: Mach number [].
% - gamma: gas adiabatic coefficient [].
% ----------------------------------------------------------------------- %
% Output:
% - A_At: corresponding A/At ratio [].
% ----------------------------------------------------------------------- %

function A_At = M2Ar(M,gamma)

A_At = sqrt(1./M.^2.*(2./(gamma+1).*(1 + (gamma-1)./2.*M.^2)).^((gamma+1)./(gamma-1)));

end

% ----------------------------------------------------------------------- %