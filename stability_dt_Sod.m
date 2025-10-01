% ----------------------------------------------------------------------- %
% Function that returns the time step required by the solver to be stable,
% according to a CFL-like condition taken from Sod, G. A. "A Survey of
% Several Finite Difference Methods for Systems of Nonlinear Hyperbolic
% Conservation Laws".
% Valid function for a 1D mesh.
% To be used with the MacCormack scheme.
% ----------------------------------------------------------------------- %
% Input:
% - u: velocity vector [m/s] [1 x N].
% - T: temperature vector [K] [1 x N].
% - Rg: gas constant vector [J/(kgK)] [1 x N].
% - gamma: gas adiabatic coefficient vector [] [1 x N].
% - sigma: CFL-like stability parameter.
% - dx: grid's cell width [m].
% ----------------------------------------------------------------------- %
% Output:
% - dt: time step required [s].
% ----------------------------------------------------------------------- %

function dt = stability_dt_Sod(u,T,Rg,gamma,sigma,dx)

% Local sound speed vector:
c = sqrt(gamma.*Rg.*T); % [m/s]

% Time step:
dt = sigma*dx/max(abs(u) + c); % [s]

end

% ----------------------------------------------------------------------- %