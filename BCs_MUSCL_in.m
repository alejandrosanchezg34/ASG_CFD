% ----------------------------------------------------------------------- %
% Function called by "BCs_MUSCL.m", and used to impose the inflow boundary
% conditions, i.e. the boundary conditions at the inlet or left boundary.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density vector [kg/m3] [1 x N+4].
% - u: velocity vector [m/s] [1 x N+4].
% - p: pressure vector [Pa] [1 x N+4].
% - Rg: gas constant [J/(kgK)] [1 x N+4].
% - gamma: gas adiabatic coefficient [] [1 x N+4].
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
% * Subscripts "in" and "out" refer to the position of the BC ghost cell.
% ----------------------------------------------------------------------- %

function [rho_L_out,u_L_out,e_L_out,T_L_out,p_L_out, ...
    rho_L_in,u_L_in,e_L_in,T_L_in,p_L_in] = BCs_MUSCL_in(rho,u,p,Rg,gamma,rho_in,u_in,p_in,BCs)

% Impose inflow boundary conditions:
switch BCs
    case {"zero_gradient"}

        % Transmissive inlet (impose zero gradient):
        rho_L_in = rho(1,3); % [kg/m3] - Density
        u_L_in = u(1,3); % [m/s] - Velocity
        p_L_in = p(1,3); % [Pa] - Pressure

        rho_L_out = rho(1,4); % [kg/m3] - Density
        u_L_out = u(1,4); % [m/s] - Velocity
        p_L_out = p(1,4); % [Pa] - Pressure

    case {"supersonic_in_free_out"}

        % Fix all variables at inlet (fix density, velocity and pressure):
        rho_L_in = 2*rho_in - rho(1,3); % [kg/m3] - Density
        u_L_in = 2*u_in - u(1,3); % [m/s] - Velocity
        p_L_in = 2*p_in - p(1,3); % [Pa] - Pressure

        % Get outer node variables (assuming transmissive conditions):
        rho_L_out = rho_L_in; % [kg/m3] - Density
        u_L_out = u_L_in; % [m/s] - Velocity
        p_L_out = p_L_in; % [Pa] - Pressure

    case {"subsonic_in_free_out","subsonic_in_fixed_out"}

        % Fix density and pressure, but allow the velocity to float:
        rho_L_in = 2*rho_in - rho(1,3); % [kg/m3] - Density
        u_L_in = 2*u(1,3) - u(1,4); % [m/s] - Velocity
        p_L_in = 2*p_in - p(1,3); % [Pa] - Pressure

        % Get outer node variables (assuming transmissive conditions):
        rho_L_out = rho_L_in; % [kg/m3] - Density
        u_L_out = u_L_in; % [m/s] - Velocity
        p_L_out = p_L_in; % [Pa] - Pressure

    otherwise
        error("Boundary conditions are not valid! Try another configuration or define new ones.")
end

% Get inlet/left approximated gas properties:
Rg_L_in = Rg(1,2); % [J/(kgK)] - Gas constant (inner halo node)
Rg_L_out = Rg(1,1); % [J/(kgK)] - Gas constant (outer halo node)
gamma_L_in = gamma(1,2); % [] - Gas adiabatic coefficient (inner halo node)
gamma_L_out = gamma(1,1); % [] - Gas adiabatic coefficient (outer halo node)
cv_L_in = getcv(Rg_L_in,gamma_L_in); % [J/(kgK)] - Gas specific heat at constant volume (inner halo node)
cv_L_out = getcv(Rg_L_out,gamma_L_out); % [J/(kgK)] - Gas specific heat at constant volume (outer halo node)

% Compute corresponding temperatures T:
T_L_in = T_from_state_equation(rho_L_in,Rg_L_in,p_L_in); % [K] - Temperature (inner halo node)
T_L_out = T_from_state_equation(rho_L_out,Rg_L_out,p_L_out); % [K] - Temperature (outer halo node)

% Compute corresponding internal specific energies e:
e_L_in = cv_L_in*T_L_in; % [m2/s2] - Specific energy (inner halo node)
e_L_out = cv_L_out*T_L_out; % [m2/s2] - Specific energy (outer halo node)

end

% ----------------------------------------------------------------------- %