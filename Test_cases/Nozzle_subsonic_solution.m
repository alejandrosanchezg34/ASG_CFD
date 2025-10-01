% ----------------------------------------------------------------------- %
% This function provides the exact analytical solution for a 1D nozzle
% flow considering it is purely subsonic, assuming a perfect gas model and
% an isentropic process along the whole nozzle.
% ----------------------------------------------------------------------- %
% Input:
% - A_f: cross-section area analytical function (function handle).
% - At: nozzle throat area [m2].
% - x: vector of grid coordinates [m] [1 x N].
% - Rg: gas constant [J/(kgK)].
% - gamma: gas adiabatic coefficient [].
% - p_0: stagnation (combustion chamber) pressure [Pa].
% - rho_0: stagnation (combustion chamber) density [kg/m3].
% - T_0: stagnation (combustion chamber) temperature [K].
% - p_out: nozzle exit pressure [Pa].
% ----------------------------------------------------------------------- %
% Output:
% - M: Mach number distribution [] [N x 1].
% - p: pressure distribution [Pa] [N x 1].
% - rho: density distribution [kg/m3] [N x 1].
% - T: temperature distribution [K] [N x 1].
% - u: velocity distribution [m/s] [N x 1].
% ----------------------------------------------------------------------- %

function [M,p,rho,T,u] = Nozzle_subsonic_solution(A_f,At,x,Rg,gamma,p_0,rho_0,T_0,p_out)

% Preliminary data:
N = length(x); % Number of grid points
A = arrayfun(A_f,x);
M_ini = 0.3; % [] - Initial guess for the Mach number
Ae = A(end); % [m2] - Exit area

% Find exit Mach number from imposed exit pressure:
fun_pe2Me = @(M) M2pr(M,gamma) - p_out/p_0; % Function to solve
Me = fzero(fun_pe2Me,M_ini,optimset('Display','off')); % [] - Exit Mach number

% Find associated sonic throat area:
Ae2Astar = M2Ar(Me,gamma); % [] - Ae/Astar ratio
Astar = Ae/Ae2Astar; % [m2] - Sonic throat area

if Astar >= At
    error("The nozzle exact solver considering purely subsonic flow has found " + ...
        "a sonic throat area higher than the actual throat area: there CAN NOT be " + ...
        "purely subsonic flow in this case. The flow is already CHOCKED.");
end

% Find Mach number for each area value:
M = zeros(N,1); % Initialize
for i = 1:N    
    fun_M2Ar = @(M) M2Ar(M,gamma) - A(i)/Astar; % Function to solve
    M(i) = fzero(fun_M2Ar,M_ini,optimset('Display','off')); % [] - Mach number
end

% Interpolate for 'NaN' values:
M(isnan(M)) = interp1(x(~isnan(M)),M(~isnan(M)),x(isnan(M)),"spline");

% Interpolate for M > 1 values:
if any(M > 1)
    M(M > 1) = interp1(x(~(M>1)),M(~(M>1)),x(M>1),"spline");
end

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