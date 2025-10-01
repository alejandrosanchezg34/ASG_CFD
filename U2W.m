% ----------------------------------------------------------------------- %
% Function that, given the solution vector U, returns the corresponding
% vector of primitive variables W.
% Valid function for 1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - U: solution vector ([3 x 1] or [3 x N]).
% - A: cross-section area [m2] (may be given as a [1 x N] vector).
% - Rg: gas constant [J/(kgK)] (may be given as a [1 x N] vector).
% - gamma: gas adiabatic coefficient [] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %
% Output:
% - W: vector of primitive variables [rho,u,p]^T [kg/m3,m/s,Pa] ([3 x 1] or
% [3 x N]).
% ----------------------------------------------------------------------- %

function W = U2W(U,A,Rg,gamma)

% Obtain flow-field variables from 'U':
[rho,u,~,~,p] = U2vars(U,A,Rg,gamma);

% Obtain corresponding vector of primitive variables 'W':
W = [rho;
     u;
     p];

end

% ----------------------------------------------------------------------- %