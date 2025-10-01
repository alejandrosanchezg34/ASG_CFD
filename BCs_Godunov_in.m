% ----------------------------------------------------------------------- %
% Function called by "BCs_Godunov.m", and used to impose the inflow
% boundary conditions, i.e. the boundary conditions at the inlet or left
% boundary.
% An exact Riemann solver or a HLLC approximate solver is used in some
% cases.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density vector [kg/m3] [1 x N+2].
% - u: velocity vector [m/s] [1 x N+2].
% - p: pressure vector [Pa] [1 x N+2].
% - Rg: gas constant [J/(kgK)] [1 x N+2].
% - gamma: gas adiabatic coefficient [] [1 x N+2].
% - rho_in: if needed, inlet density to be fixed [kg/m3].
% - u_in: if needed, inlet velocity to be fixed [m/s].
% - p_in: if needed, inlet pressure to be fixed [Pa].
% - BCs: string specifying the type of boundary conditions to impose.
% - Riemann_solver: Riemann problem solver to be used ("HLLC" or "exact").
% ----------------------------------------------------------------------- %
% Output:
% - rho_L: inlet/left density [kg/m3].
% - u_L: inlet/left velocity [m/s].
% - e_L: inlet/left specific energy [m2/s2].
% - T_L: inlet/left temperature [K].
% - p_L: inlet/left pressure [Pa].
% ----------------------------------------------------------------------- %

function [rho_L,u_L,e_L,T_L,p_L] = BCs_Godunov_in(rho,u,p,Rg,gamma,rho_in,u_in,p_in,BCs,Riemann_solver)

% Impose inflow boundary conditions:
switch BCs
    case {"zero_gradient"}

        % Transmissive inlet (impose zero gradient):
        rho_L = rho(1,2); % [kg/m3] - Density
        u_L = u(1,2); % [m/s] - Velocity
        p_L = p(1,2); % [Pa] - Pressure

    case {"supersonic_in_free_out"}

        % Imposing a Riemann solver at the inlet to fix all variables
        % (density, velocity and pressure):
        
        % System of equations to solve:
        % z(1): rho
        % z(2): u
        % z(3): p
        
        % Prepare left and right data for the Riemann solver:
        WL = @(z) [z(1); z(2); z(3)]; % Left (node 1 data is unknown)
        WR = [rho(1,2); u(1,2); p(1,2)]; % Right (node 2 data is known)

        % Define intercell gas properties:
        Rg_int = harmmean(Rg(1,1:2)); % [J/(kgK)] - Gas constant
        gamma_int = harmmean(gamma(1,1:2)); % [] - Gas adiabatic coefficient

        % Choose Riemann solver to solve system of equations (and to meet BCs):
        switch Riemann_solver
            case "HLLC"
                W_BCs = @(z) HLLC2WSol(WL(z),WR,Rg_int,gamma_int,BCs) - [rho_in; u_in; p_in]; % Desired data at x = 0 (left boundary)
            case "exact"
                W_BCs = @(z) exact2WSol(WL(z),WR,gamma_int,BCs) - [rho_in; u_in; p_in]; % Desired data at x = 0 (left boundary)
        end

        z = fsolve(W_BCs,[rho_in; u_in; p_in],optimoptions('fsolve','Display','off')); % Find corresponding left data

        % Assign variables:
        rho_L = z(1); % [kg/m3] - Density
        u_L = z(2); % [m/s] - Velocity
        p_L = z(3); % [Pa] - Pressure

    case {"subsonic_in_free_out","subsonic_in_fixed_out"}

        % Imposing a Riemann solver at the inlet to fix density and
        % pressure (and allow the velocity to float):
        
        % System of equations to solve:
        % z(1): rho
        % z(2): u
        % z(3): p
        
        % Prepare left and right data for the Riemann solver:
        WL = @(z) [z(1); z(2); z(3)]; % Left (node 1 data is unknown)
        WR = [rho(1,2); u(1,2); p(1,2)]; % Right (node 2 data is known)

        % Define intercell gas properties:
        Rg_int = harmmean(Rg(1,1:2)); % [J/(kgK)] - Gas constant
        gamma_int = harmmean(gamma(1,1:2)); % [] - Gas adiabatic coefficient

        % Choose Riemann solver to solve system of equations (and to meet BCs):
        switch Riemann_solver
            case "HLLC"
                W_BCs = @(z) HLLC2WSol(WL(z),WR,Rg_int,gamma_int,BCs) - [rho_in; p_in]; % Desired data at x = 0 (left boundary)
            case "exact"
                W_BCs = @(z) exact2WSol(WL(z),WR,gamma_int,BCs) - [rho_in; p_in]; % Desired data at x = 0 (left boundary)
        end

        z = fsolve(W_BCs,[rho_in; u(1,2); p_in],optimoptions('fsolve','Display','off', ...
                    'Algorithm','levenberg-marquardt')); % Find corresponding left data

        % Assign variables:
        rho_L = z(1); % [kg/m3] - Density
        u_L = z(2); % [m/s] - Velocity
        p_L = z(3); % [Pa] - Pressure

    otherwise
        error("Boundary conditions are not valid! Try another configuration or define new ones.")
end

% Get inlet/left approximated gas properties:
Rg_L = Rg(1,1); % [J/(kgK)] - Gas constant
gamma_L = gamma(1,1); % [] - Gas adiabatic coefficient
cv_L = getcv(Rg_L,gamma_L); % [J/(kgK)] - Gas specific heat at constant volume

% Compute corresponding temperature T:
T_L = T_from_state_equation(rho_L,Rg_L,p_L); % [K] - Temperature

% Compute corresponding internal specific energy e:
e_L = cv_L*T_L; % [m2/s2] - Specific energy

end

% ----------------------------------------------------------------------- %

%% Other functions

% ----------------------------------------------------------------------- %
% Function required to obtain WSol from the HLLC Riemann solver:
function WSol = HLLC2WSol(WL,WR,Rg,gamma,BCs)

% Avoid negative densities and pressures and complex states by penalizing these cases:
if (any(WL([1 3],1) < 0) || any(WR([1 3],1) < 0) || ~isreal(WL) || ~isreal(WR))
    Uint = NaN(3,1);
else
    % Obtain solution vector Uint:
    [~,~,~,Uint] = HLLC_solver(WL,WR,Rg,gamma);
end

% Obtain corresponding vector of primitive variables:
switch BCs
    case {"supersonic_in_free_out"}
        WSol = U2W(Uint,1,Rg,gamma); % [rho; u; p]

    case {"subsonic_in_free_out","subsonic_in_fixed_out"}
        W = U2W(Uint,1,Rg,gamma); % [rho; u; p]
        WSol = [W(1); W(3)]; % [rho; p]
end

end
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Function required to obtain WSol from the exact Riemann solver:
function WSol = exact2WSol(WL,WR,gamma,BCs)

% Avoid negative densities and pressures and complex states by penalizing these cases:
if (any(WL([1 3],1) < 0) || any(WR([1 3],1) < 0) || ~isreal(WL) || ~isreal(WR))
    W = NaN(3,1);
else
    % Obtain original WSol vector:
    W = Riemann_solution(0,1,gamma,WL,WR);
end

% Obtain corresponding vector of primitive variables:
switch BCs
    case {"supersonic_in_free_out"}
        WSol = W; % [rho; u; p]

    case {"subsonic_in_free_out","subsonic_in_fixed_out"}
        WSol = [W(1); W(3)]; % [rho; p]
end

end
% ----------------------------------------------------------------------- %