% ----------------------------------------------------------------------- %
% Function that, given the solution vector U, returns the corresponding
% flux vector F.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - U: solution vector ([3 x 1] or [3 x N]).
% - A: cross-section area [m2] (may be given as a [1 x N] vector).
% - Rg: gas constant [J/(kgK)] (may be given as a [1 x N] vector).
% - gamma: gas adiabatic coefficient [] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %
% Output:
% - F: flux vector ([3 x 1] or [3 x N]).
% ----------------------------------------------------------------------- %

function F = U2F(U,A,Rg,gamma)

% Get flow-field variables from 'U':
[rho,u,e,~,p] = U2vars(U,A,Rg,gamma);

% Get corresponding flux vector 'F':
F = vars2F(rho,u,e,p,A);

end

% ----------------------------------------------------------------------- %