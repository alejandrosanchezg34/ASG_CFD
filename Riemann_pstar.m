% ----------------------------------------------------------------------- %
% Function that implements the Newton-Raphson iterative method to find the
% root of Riemann's algebraic equation, i.e. the pressure p* at the Star
% Region.
% ----------------------------------------------------------------------- %
% Input:
% - gamma: gas ratio of specific heats.
% - WL: vector of left primitive variables [rhoL,uL,pL]^T [kg/m3,m/s,Pa].
% - WR: vector of right primitive variables [rhoR,uR,pR]^T [kg/m3,m/s,Pa].
% ----------------------------------------------------------------------- %
% Output:
% - pstar: pressure solution, i.e. pressure at the Star Region [Pa].
% - it: number of iterations carried out.
% ----------------------------------------------------------------------- %

function [pstar,it] = Riemann_pstar(gamma,WL,WR)

% Getting left and right pressures:
pL = WL(3); % [Pa]
pR = WR(3); % [Pa]

% Numerical data:
tol = 1e-6; % Solver tolerance
Nmax = 1e2; % Number of maximum iterations
p_old = 0.5*(pL+pR); % [Pa] - Initial guess

% Initializing solver:
it = 0; % Number of iterations
cvg = 0; % Boolean to indicate whether convergence has been achieved or not

% Newton-Raphson method:
while it < Nmax
    it = it + 1; % New iteration

    % Obtain f(p_old) and f'(p_old):
    [f,fp] = Riemann_f_p(gamma,WL,WR,p_old);

    % New pressure value:
    p_new = p_old - f/fp; % [Pa]

    % Relative pressure change:
    CHA = abs(p_new - p_old)/(0.5*(p_new + p_old));

    % Check convergence:
    if CHA < tol
        cvg = 1; % Converged!
        break
    end

    % Next iteration:
    p_old = p_new;
end

% Return final pressure value:
pstar = p_new; % [Pa]

if pstar < 0
    disp("Negative p* has been found.");
end

if ~cvg
    disp("Riemann p* solver did not converge! p* has not been found.");
end

end

% ----------------------------------------------------------------------- %