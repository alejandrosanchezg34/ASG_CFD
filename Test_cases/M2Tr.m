% ----------------------------------------------------------------------- %
% This function provides the temperature ratio T/T0 associated to a given
% Mach number M, where T0 is the stagnation temperature.
% ----------------------------------------------------------------------- %
% Input:
% - M: Mach number [].
% - gamma: gas adiabatic coefficient [].
% ----------------------------------------------------------------------- %
% Output:
% - T_T0: corresponding T/T0 ratio [].
% ----------------------------------------------------------------------- %

function T_T0 = M2Tr(M,gamma)

T_T0 = (1 + (gamma-1)./2.*M.^2).^(-1);

end

% ----------------------------------------------------------------------- %