% ----------------------------------------------------------------------- %
% Function called by "BCs_MacCormack.m", and used to impose the inflow
% boundary conditions, i.e. the boundary conditions at the inlet or left
% boundary.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density vector [kg/m3] [1 x N].
% - u: velocity vector [m/s] [1 x N].
% - p: pressure vector [Pa] [1 x N].
% - Rg: gas constant [J/(kgK)] [1 x N].
% - gamma: gas adiabatic coefficient [] [1 x N].
% - rho_in: if needed, inlet density to be fixed [kg/m3].
% - u_in: if needed, inlet velocity to be fixed [m/s].
% - p_in: if needed, inlet pressure to be fixed [Pa].
% - BCs: string specifying the type of boundary conditions to impose.
% ----------------------------------------------------------------------- %
% Output:
% - rho_L: inlet/left density [kg/m3].
% - u_L: inlet/left velocity [m/s].
% - e_L: inlet/left specific energy [m2/s2].
% - T_L: inlet/left temperature [K].
% - p_L: inlet/left pressure [Pa].
% ----------------------------------------------------------------------- %

function [rho_L,u_L,e_L,T_L,p_L] = BCs_MacCormack_in(rho,u,p,Rg,gamma,rho_in,u_in,p_in,BCs)

% Impose inflow boundary conditions:
switch BCs
    case {"zero_gradient"}

        % Transmissive inlet (impose zero gradient):
        rho_L = rho(1,2); % [kg/m3] - Density
        u_L = u(1,2); % [m/s] - Velocity
        p_L = p(1,2); % [Pa] - Pressure

    case {"supersonic_in_free_out"}

        % Fix all variables at inlet (fix density, velocity and pressure):
        rho_L = rho_in; % [kg/m3] - Density
        u_L = u_in; % [m/s] - Velocity
        p_L = p_in; % [Pa] - Pressure

    case {"subsonic_in_free_out","subsonic_in_fixed_out"}
        
        % Fix density and pressure, but allow the velocity to float:
        rho_L = rho_in; % [kg/m3] - Density
        u_L = 2*u(1,2) - u(1,3); % [m/s] - Velocity
        p_L = p_in; % [Pa] - Pressure

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