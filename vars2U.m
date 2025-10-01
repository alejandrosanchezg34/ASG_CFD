% ----------------------------------------------------------------------- %
% Function that, given the flow-field variables, returns the corresponding
% solution vector U.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density [kg/m3].
% - u: velocity [m/s].
% - e: internal specific energy [m2/s2].
% - A: cross-section area [m2].
% * May be given as vectors [1 x N].
% ----------------------------------------------------------------------- %
% Output:
% - U: solution vector ([3 x 1] or [3 x N]).
% ----------------------------------------------------------------------- %

function U = vars2U(rho,u,e,A)

% Arrange U vector:
U = [rho.*A;
     rho.*A.*u;
     rho.*(e + 0.5*u.^2).*A];

end

% ----------------------------------------------------------------------- %