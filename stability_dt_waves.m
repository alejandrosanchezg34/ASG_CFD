% ----------------------------------------------------------------------- %
% Function that returns the time step required by the solver to be stable,
% according to a CFL-like condition taken from Toro, E. F. "Riemann Solvers
% and Numerical Methods for Fluid Dynamics: A Practical Introduction".
% Valid function for a 1D grid and for the Godunov method.
% ----------------------------------------------------------------------- %
% Input:
% - SL: vector of speeds associated to the left non-linear waves [m/s]
% (from a 'i+1/2' Riemann problem) [N+1 x 1].
% - SR: vector of speeds associated to the right non-linear waves [m/s]
% (from a 'i+1/2' Riemann problem) [N+1 x 1].
% - dx: grid's cell width [m].
% - sigma: stability parameter.
% ----------------------------------------------------------------------- %
% Output:
% - dt: time step required [s].
% ----------------------------------------------------------------------- %

function dt = stability_dt_waves(SL,SR,dx,sigma)

% Initializing:
Si_max = zeros(length(SL),1);

% Maximum non-linear wave speed present in the solution of the Riemann
% problem centered at 'i+1/2':
for i = 1:length(SL)
    Si_max(i) = max(abs(SL(i)),abs(SR(i)));
end

% Maximum non-linear wave speed throughout the 1D domain:
Sn_max = max(Si_max);

% Time step required:
dt = sigma*dx/Sn_max; % [s]

end

% ----------------------------------------------------------------------- %