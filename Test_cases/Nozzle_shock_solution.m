% ----------------------------------------------------------------------- %
% This function provides the exact analytical solution for a 1D nozzle
% flow, assuming a perfect gas model and an isentropic process along the
% whole nozzle, but considering there is a shock wave inside.
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
% - p_out: nozzle exit pressure [Pa].
% ----------------------------------------------------------------------- %
% Output:
% - M: Mach number distribution [] [N x 1].
% - p: pressure distribution [Pa] [N x 1].
% - rho: density distribution [kg/m3] [N x 1].
% - T: temperature distribution [K] [N x 1].
% - u: velocity distribution [m/s] [N x 1].
% - p_0_vec: stagnation pressures vector [Pa] [N x 1].
% - rho_0_vec: stagnation densities vector [kg/m3] [N x 1].
% - T_0_vec: stagnation temperatures vector [K] [N x 1].
% ----------------------------------------------------------------------- %

function [M,p,rho,T,u,p_0_vec,rho_0_vec,T_0_vec] = Nozzle_shock_solution(A_f,dAdx_f,At,x,Rg,gamma,p_0,rho_0,T_0,p_out)

% Number of grid points:
N = length(x);

% Build 'A' and 'dAdx' vectors:
A = arrayfun(A_f,x);
dAdx = arrayfun(dAdx_f,x);

% Initializing matrices:
M = zeros(N,1);
p_0_vec = zeros(N,1);
rho_0_vec = zeros(N,1);
T_0_vec = zeros(N,1);

% Initial guesses:
M_sub = 0.5; % [] - Subsonic Mach number
M_sup = 3.0; % [] - Supersonic Mach number

% pe*Ae/(p0e*Ate) ratio:
peAe_p0eAte = p_out/p_0*A(end)/At; % []

% Find exit Mach number due to shock:
fun_Me = @(M) M2pr(M,gamma)*M2Ar(M,gamma) - peAe_p0eAte;
Me = fzero(fun_Me,M_sub,optimset('Display','off')); % []

% Find exit pressure ratio:
pe_p0e = M2pr(Me,gamma); % []

% Stagnation pressure ratio across normal shock:
p02_p01 = p_out/p_0/pe_p0e; % []

% Mach number in front of shock:
fun_M1 = @(M) ((((gamma+1)*M^2)/((gamma-1)*M^2 + 2))^(gamma/(gamma-1))*((gamma+1)/(2*gamma*M^2 - ...
    (gamma-1)))^(1/(gamma-1))) - p02_p01;
M1 = fzero(fun_M1,M_sup,optimset('Display','off')); % []

% Find area ratio where the shock wave is located:
A1_At = M2Ar(M1,gamma); % []

if (A1_At > A(end)/At)
    error("For the given exit pressure, there is no shock wave INSIDE the nozzle.");
end

% Find area ratio at exit and after-shock sonic throat area:
Ae_Astar = M2Ar(Me,gamma); % [] - Exit area ratio
At_after = A(end)/Ae_Astar; % [m2] - After-shock sonic throat area

% ---
% Continue by computing the analytical solution over A(x):
% ---

% Boolean defining whether the shock wave has been found or not:
shock = 0;

% Find Mach number for each area value:
M_ini = M_sub; % Initial guess for convergent region

for i = 1:N
    Astar = At; % Sonic throat area [m2]

    if (dAdx(i) >= 0 && A(i)/At < A1_At)
        M_ini = M_sup; % Initial guess for divergent region BEFORE shock

    elseif (dAdx(i) >= 0 && A(i)/At >= A1_At)
        M_ini = M_sub; % Initial guess for divergent region AFTER shock
        Astar = At_after; % Sonic throat area [m2]
        
        % In case the shock is found for the first time:
        if ~shock
            shock = 1; % Shock has been found
            ind_shock = i; % Index associated to shock location
        end
    end

    fun_M2A = @(M) M2Ar(M,gamma) - A(i)/Astar; % Function to solve
    M(i) = fzero(fun_M2A,M_ini,optimset('Display','off')); % [] - Mach number
end

% Interpolate for 'NaN' values:
M(isnan(M)) = interp1(x(~isnan(M)),M(~isnan(M)),x(isnan(M)),"spline");

% Interpolate for M < 0 values:
if any(M < 0)
    M(M < 0) = interp1(x(~(M<0)),M(~(M<0)),x(M<0),"spline");
end

% After-shock stagnation variables:
T_0_2 = T_0; % [K]
p_0_2 = p_0*p02_p01; % [Pa]
rho_0_2 = p_0_2/(Rg*T_0_2); % [kg/m3]

% Stagnation variable vectors:
for i = 1:N
    if i < ind_shock % Before shock
        T_0_vec(i) = T_0; % [K]
        p_0_vec(i) = p_0; % [Pa]
        rho_0_vec(i) = rho_0; % [kg/m3]
    
    else % After shock
        T_0_vec(i) = T_0_2; % [K]
        p_0_vec(i) = p_0_2; % [Pa]
        rho_0_vec(i) = rho_0_2; % [kg/m3]
    end
end

% Find the distribution of thermodynamic variables:
p = M2pr(M,gamma).*p_0_vec; % [Pa] - Pressure
rho = M2rhor(M,gamma).*rho_0_vec; % [kg/m3] - Density
T = M2Tr(M,gamma).*T_0_vec; % [K] - Temperature
a = sqrt(gamma*Rg.*T); % [m/s] - Sound speed
u = M.*a; % [m/s] - Flow velocity

end

% ----------------------------------------------------------------------- %