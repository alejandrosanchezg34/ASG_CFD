% ----------------------------------------------------------------------- %
% Function called by "BCs_Godunov.m", and used to impose the outflow
% boundary conditions, i.e. the boundary conditions at the outlet or right
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
% - p_out: if needed, outlet pressure to be fixed [Pa].
% - BCs: string specifying the type of boundary conditions to impose.
% - Riemann_solver: Riemann problem solver to be used ("HLLC" or "exact").
% ----------------------------------------------------------------------- %
% Output:
% - rho_R: outlet/right density [kg/m3].
% - u_R: outlet/right velocity [m/s].
% - e_R: outlet/right specific energy [m2/s2].
% - T_R: outlet/right temperature [K].
% - p_R: outlet/right pressure [Pa].
% ----------------------------------------------------------------------- %

function [rho_R,u_R,e_R,T_R,p_R] = BCs_Godunov_out(rho,u,p,Rg,gamma,p_out,BCs,Riemann_solver)

% Impose outflow boundary conditions:
switch BCs
    case {"zero_gradient","supersonic_in_free_out","subsonic_in_free_out"}

        % Transmissive outlet (impose zero gradient):
        rho_R = rho(1,end-1); % [kg/m3] - Density
        u_R = u(1,end-1); % [m/s] - Velocity
        p_R = p(1,end-1); % [Pa] - Pressure

    case {"subsonic_in_fixed_out"}

        % Imposing a Riemann solver at the outlet to fix only pressure (and
        % allowing density and velocity to float):
        
        % System of equations to solve:
        % z(1): rho
        % z(2): u
        % z(3): p
        
        % Prepare left and right data for the Riemann solver:
        WL = [rho(1,end-1); u(1,end-1); p(1,end-1)]; % Left (node N+1 data is known)
        WR = @(z) [z(1); z(2); z(3)]; % Right (node N+2 data is unknown)

        % Define intercell gas properties:
        Rg_int = harmmean(Rg(1,end-1:end)); % [J/(kgK)] - Gas constant
        gamma_int = harmmean(gamma(1,end-1:end)); % [] - Gas adiabatic coefficient

        % Choose Riemann solver to solve system of equations (and to meet BCs):
        switch Riemann_solver
            case "HLLC"
                W_BCs = @(z) HLLC2WSol(WL,WR(z),Rg_int,gamma_int,BCs) - p_out; % Desired data at x = L (right boundary)
            case "exact"
                W_BCs = @(z) exact2WSol(WL,WR(z),gamma_int,BCs) - p_out; % Desired data at x = L (right boundary)
        end

        z = fsolve(W_BCs,[rho(1,end-1); u(1,end-1); p_out],optimoptions('fsolve','Display','off', ...
                    'Algorithm','levenberg-marquardt')); % Find corresponding right data

        % Assign variables:
        rho_R = z(1); % [kg/m3] - Density
        u_R = z(2); % [m/s] - Velocity
        p_R = z(3); % [Pa] - Pressure

    otherwise
        error("Boundary conditions are not valid! Try another configuration or define new ones.")
end

% Get outlet/right approximated gas properties:
Rg_R = Rg(1,end); % [J/(kgK)] - Gas constant
gamma_R = gamma(1,end); % [] - Gas adiabatic coefficient
cv_R = getcv(Rg_R,gamma_R); % [J/(kgK)] - Gas specific heat at constant volume

% Compute corresponding temperature T:
T_R = T_from_state_equation(rho_R,Rg_R,p_R); % [K] - Temperature

% Compute corresponding internal specific energy e:
e_R = cv_R*T_R; % [m2/s2] - Specific energy

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
    case {"subsonic_in_fixed_out"}
        W = U2W(Uint,1,Rg,gamma); % [rho; u; p]
        WSol = W(3); % [p]
end

end
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Function required to obtain WSol from the exact Riemann solver:
function WSol = exact2WSol(WL,WR,gamma,BCs)

% Avoid negative and complex states by penalizing these cases:
if (any(WL([1 3],1) < 0) || any(WR([1 3],1) < 0) || ~isreal(WL) || ~isreal(WR))
    W = NaN(3,1);
else
    % Obtain original WSol vector:
    W = Riemann_solution(0,1,gamma,WL,WR);
end

% Obtain corresponding vector of primitive variables:
switch BCs
    case {"subsonic_in_fixed_out"}
        WSol = W(3); % [p]
end

end
% ----------------------------------------------------------------------- %