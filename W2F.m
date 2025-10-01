% ----------------------------------------------------------------------- %
% Function that, given the Riemann problem solution vector W, returns the
% corresponding flux vector F.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - W: solution vector ([3 x 1] or [3 x N]).
% - A: cross-section area [m2] (may be given as a [1 x N] vector).
% - Rg: gas constant [J/(kgK)] (may be given as a [1 x N] vector).
% - gamma: gas adiabatic coefficient [] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %
% Output:
% - F: flux vector ([3 x 1] or [3 x N]).
% ----------------------------------------------------------------------- %

function F = W2F(W,A,Rg,gamma)

% Get flow-field variables from 'W':
[rho,u,e,~,p] = W2vars(W,Rg,gamma);

% Get corresponding flux vector 'F':
F = vars2F(rho,u,e,p,A);

end

% ----------------------------------------------------------------------- %