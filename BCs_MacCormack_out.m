% ----------------------------------------------------------------------- %
% Function called by "BCs_MacCormack.m", and used to impose the outflow
% boundary conditions, i.e. the boundary conditions at the outlet or right
% boundary.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density vector [kg/m3] [1 x N].
% - u: velocity vector [m/s] [1 x N].
% - p: pressure vector [Pa] [1 x N].
% - Rg: gas constant [J/(kgK)] [1 x N].
% - gamma: gas adiabatic coefficient [] [1 x N].
% - p_out: if needed, outlet pressure to be fixed [Pa].
% - BCs: string specifying the type of boundary conditions to impose.
% ----------------------------------------------------------------------- %
% Output:
% - rho_R: outlet/right density [kg/m3].
% - u_R: outlet/right velocity [m/s].
% - e_R: outlet/right specific energy [m2/s2].
% - T_R: outlet/right temperature [K].
% - p_R: outlet/right pressure [Pa].
% ----------------------------------------------------------------------- %

function [rho_R,u_R,e_R,T_R,p_R] = BCs_MacCormack_out(rho,u,p,Rg,gamma,p_out,BCs)

% Impose outflow boundary conditions:
switch BCs
    case {"zero_gradient"}

        % Transmissive outlet (impose zero gradient):
        rho_R = rho(1,end-1); % [kg/m3] - Density
        u_R = u(1,end-1); % [m/s] - Velocity
        p_R = p(1,end-1); % [Pa] - Pressure

    case {"supersonic_in_free_out","subsonic_in_free_out"}

        % Unconstrained outlet (let all variables float):
        rho_R = 2*rho(1,end-1) - rho(1,end-2); % [kg/m3] - Density
        u_R = 2*u(1,end-1) - u(1,end-2); % [m/s] - Velocity
        p_R = 2*p(1,end-1) - p(1,end-2); % [Pa] - Pressure

    case "subsonic_in_fixed_out"

        % Fix outlet pressure, but let the rest of variables float:
        rho_R = 2*rho(1,end-1) - rho(1,end-2); % [kg/m3] - Density
        u_R = 2*u(1,end-1) - u(1,end-2); % [m/s] - Velocity
        p_R = p_out; % [Pa] - Pressure

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