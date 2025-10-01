% ----------------------------------------------------------------------- %
% This function provides the exact analytic solution for a Prandtl-Meyer
% expansion wave, given the free stream flow variables as well as the turn
% angle.
% ----------------------------------------------------------------------- %
% Input:
% - p1: free stream pressure [Pa].
% - T1: free stream temperature [K].
% - M1: free stream Mach number [].
% - Rg: gas constant [J/(kgK)].
% - gamma: gas adiabatic coefficient [].
% - theta: turn angle [deg].
% - N: number of steps to solve the problem (if nothing is given, 1 value
% is taken, which is for the given theta angle).
% ----------------------------------------------------------------------- %
% Output:
% - M2: after-expansion Mach number [].
% - p2: after-expansion pressure [Pa].
% - rho2: after-expansion density [kg/m3].
% - T2: after-expansion temperature [K].
% - u2: after-expansion velocity [m/s].
% - ang: vector of evaluated turn angles [deg].
% - it: vector of required number of iterations.
% All these variables are given as vectors over a continuous and linear
% variation of the angle theta, starting at 0 deg (free stream).
% ----------------------------------------------------------------------- %

function [M2,p2,rho2,T2,u2,ang,it] = Prandtl_Meyer_solution(p1,T1,M1,Rg,gamma,theta,N)

Nit = 1e2; % Maximum number of iterations
tol = 1e-8; % Tolerance

if nargin == 6
    N = 1; % Solve only for given theta angle
end

% Convert theta from deg to rad:
theta = deg2rad(theta); % [rad]

% Find nu(M1):
nu_M1 = PM_func(M1,gamma); % [rad]

% --- Compute maximum turning angle:
nu_max = pi/2*(sqrt((gamma+1)/(gamma-1)) - 1); % [rad]
theta_max = nu_max - nu_M1; % [rad] - Maximum turning angle

if theta > theta_max
    disp("Given turning angle for the Prandtl-Meyer expansion is higher than the maximum one!");
    p2 = NaN; T2 = NaN; M2 = NaN; rho2 = NaN; u2 = NaN; ang = NaN; it = NaN;
    return;
end
% ---

% Solve expansion for different turn angles from 0 to given theta:
ang = linspace(0,theta,N); % [rad]

M2 = zeros(N,1);
T2 = zeros(N,1);
p2 = zeros(N,1);
rho2 = zeros(N,1);
u2 = zeros(N,1);
it = zeros(N,1);

for i = 1:N
    % Find nu(M2):
    nu_M2 = ang(i) + nu_M1; % [rad]

    % Find M2 (solve Prandtl-Meyer function by bisection):
    % ------------------------------------------------------------------- %
    % Bisection limits:
    M_low = 1;
    M_up = 25;

    % First guess:
    it(i) = 1; % Iteration number
    M = (M_low + M_up)/2; % Mach number

    % Bisection solver:
    while it(i) <= Nit
        
        % Find nu(M):
        nu_M = PM_func(M,gamma);

        % If converged:
        if abs(nu_M - nu_M2) < tol
            M2(i) = M; % M2 found!
            break;
        else
            if nu_M - nu_M2 > 0
                M_up = M;
            else
                M_low = M;
            end
            it(i) = it(i) + 1; % Next iteration
            M = (M_low + M_up)/2; % Next iteration Mach number
        end
    end

    % If not converged:
    if it(i) > Nit
        disp("The bisection method to solve the Prandtl-Meyer function did not converge!");
        p2(i) = NaN; T2(i) = NaN; M2(i) = NaN; rho2(i) = NaN; u2(i) = NaN; ang(i) = NaN;
        break;
    end

    % ------------------------------------------------------------------- %

    % Find after-expansion flow variables:
    T_ratio = (1 + (gamma-1)/2*M2(i)^2)/(1 + (gamma-1)/2*M1^2);
    T2(i) = T1/T_ratio; % [K]
    p2(i) = p1/T_ratio^(gamma/(gamma-1)); % [Pa]
    rho2(i) = p2(i)/(Rg*T2(i)); % [kg/m3]
    u2(i) = M2(i)*sqrt(gamma*Rg*T2(i)); % [m/s]
end

% Return angles vector in degrees:
ang = rad2deg(ang)'; % [deg]

end

% ----------------------------------------------------------------------- %