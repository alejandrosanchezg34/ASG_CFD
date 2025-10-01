% ----------------------------------------------------------------------- %
% Function that, given the flow-field variables, returns the corresponding
% flux vector F.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density [kg/m3].
% - u: velocity [m/s].
% - e: internal specific energy [m2/s2].
% - p: pressure [Pa].
% - A: cross-section area [m2].
% * May be given as vectors [1 x N].
% ----------------------------------------------------------------------- %
% Output:
% - F: flux vector ([3 x 1] or [3 x N]).
% ----------------------------------------------------------------------- %

function F = vars2F(rho,u,e,p,A)

% Arrange F vector:
F = [rho.*A.*u;
     rho.*A.*u.^2 + p.*A;
     rho.*(e + 0.5*u.^2).*A.*u + p.*A.*u];

end

% ----------------------------------------------------------------------- %