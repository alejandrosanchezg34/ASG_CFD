% ----------------------------------------------------------------------- %
% Function that checks the convergence of the numerical solution in order
% to assess whether the steady state has been reached or not. This applies
% to transient algorithms.
% The function does so by computing the residual of the density.
% ----------------------------------------------------------------------- %
% Input:
% - rho: complete density field [kg/m3] (Nt x Nx).
% * Nt: number of time steps.
% * Nx: number of spatial steps (cells/nodes).
% ----------------------------------------------------------------------- %
% Output:
% - flag: boolean variable that states whether convergence has been met or
% not.
% ----------------------------------------------------------------------- %

function flag = check_convergence(rho)

% Set residual tolerance:
TOL = 1E-6; % [kg/m3]

% Compute density residuals:
res = abs(rho(end,:) - rho(end-1,:)); % [kg/m3]

% Check if the maximum residual is below the tolerance:
if max(res) < TOL
    flag = true;
else
    flag = false;
end

end

% ----------------------------------------------------------------------- %