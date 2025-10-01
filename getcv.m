% ----------------------------------------------------------------------- %
% Function that, given 'Rg' and 'gamma', returns 'cv'.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - Rg: gas constant [J/(kgK)] (may be given as a [1 x N] vector).
% - gamma: gas adiabatic coefficient [] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %
% Output:
% - cv: gas specific heat at constant volume [J/(kgK)] (may be given as a
% [1 x N] vector).
% ----------------------------------------------------------------------- %

function cv = getcv(Rg,gamma)

% Get 'cv' from the following relation:
cv = Rg./(gamma - 1);

end

% ----------------------------------------------------------------------- %