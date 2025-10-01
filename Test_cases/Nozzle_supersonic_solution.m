% ----------------------------------------------------------------------- %
% This function provides the exact analytical solution for a 1D nozzle
% flow, assuming a perfect gas model and an isentropic process along the
% whole nozzle, i.e. the isentropic solution is given.
% ----------------------------------------------------------------------- %
% Input:
% - A_f: cross-section area analytical function (function handle).
% - dAdx_f: cross-section area derivative analytical function (function handle).
% - At: nozzle throat area [m2].
% - x: vector of grid coordinates [m] [1 x N].
% - Rg: gas constant [J/(kgK)].
% - gamma: gas adiabatic coefficient [].
% - p_0: stagnation (combustion chamber) pressure [Pa].
% - rho_0: stagnation (combustion chamber) density [kg/m3].
% - T_0: stagnation (combustion chamber) temperature [K].
% ----------------------------------------------------------------------- %
% Output:
% - M: Mach number distribution [] [N x 1].
% - p: pressure distribution [Pa] [N x 1].
% - rho: density distribution [kg/m3] [N x 1].
% - T: temperature distribution [K] [N x 1].
% - u: velocity distribution [m/s] [N x 1].
% ----------------------------------------------------------------------- %

function [M,p,rho,T,u] = Nozzle_supersonic_solution(A_f,dAdx_f,At,x,Rg,gamma,p_0,rho_0,T_0)

% Number of grid points:
N = length(x);

% Build 'A' and 'dAdx' vectors:
A = arrayfun(A_f,x);
dAdx = arrayfun(dAdx_f,x);

% Initial guesses:
M_sub = 0.5; % [] - Subsonic Mach number
M_sup = 3.0; % [] - Supersonic Mach number

% Find Mach number for each area value:
M = zeros(N,1); % Initialize
M_ini = M_sub; % Initial guess for convergent region

for i = 1:N    
    fun_M2A = @(M) M2Ar(M,gamma) - A(i)/At; % Function to solve
    
    if (dAdx(i) >= 0)
        M_ini = M_sup; % Initial guess for divergent region
    end

    M(i) = fzero(fun_M2A,M_ini,optimset('Display','off')); % [] - Mach number
end

% Interpolate for 'NaN' values:
M(isnan(M)) = interp1(x(~isnan(M)),M(~isnan(M)),x(isnan(M)),"spline");

% Interpolate for M < 0 values:
if any(M < 0)
    M(M < 0) = interp1(x(~(M<0)),M(~(M<0)),x(M<0),"spline");
end

% Find the distribution of thermodynamic variables:
p = M2pr(M,gamma)*p_0; % [Pa] - Pressure
rho = M2rhor(M,gamma)*rho_0; % [kg/m3] - Density
T = M2Tr(M,gamma)*T_0; % [K] - Temperature
a = sqrt(gamma*Rg.*T); % [m/s] - Sound speed
u = M.*a; % [m/s] - Flow velocity

end

% ----------------------------------------------------------------------- %