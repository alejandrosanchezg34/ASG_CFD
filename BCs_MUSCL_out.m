% ----------------------------------------------------------------------- %
% Function called by "BCs_MUSCL.m", and used to impose the outflow boundary
% conditions, i.e. the boundary conditions at the outlet or right boundary.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density vector [kg/m3] [1 x N+4].
% - u: velocity vector [m/s] [1 x N+4].
% - p: pressure vector [Pa] [1 x N+4].
% - Rg: gas constant [J/(kgK)] [1 x N+4].
% - gamma: gas adiabatic coefficient [] [1 x N+4].
% - p_out: if needed, outlet pressure to be fixed [Pa].
% - BCs: string specifying the type of boundary conditions to impose.
% ----------------------------------------------------------------------- %
% Output:
% - rho_R: outlet/right density [kg/m3].
% - u_R: outlet/right velocity [m/s].
% - e_R: outlet/right specific energy [m2/s2].
% - T_R: outlet/right temperature [K].
% - p_R: outlet/right pressure [Pa].
% * Subscripts "in" and "out" refer to the position of the BC ghost cell.
% ----------------------------------------------------------------------- %

function [rho_R_in,u_R_in,e_R_in,T_R_in,p_R_in, ...
    rho_R_out,u_R_out,e_R_out,T_R_out,p_R_out] = BCs_MUSCL_out(rho,u,p,Rg,gamma,p_out,BCs)

% Impose outflow boundary conditions:
switch BCs
    case {"zero_gradient","supersonic_in_free_out","subsonic_in_free_out"}

        % Transmissive outlet (impose zero gradient):
        rho_R_in = rho(1,end-2); % [kg/m3] - Density
        u_R_in = u(1,end-2); % [m/s] - Velocity
        p_R_in = p(1,end-2); % [Pa] - Pressure

        rho_R_out = rho(1,end-3); % [kg/m3] - Density
        u_R_out = u(1,end-3); % [m/s] - Velocity
        p_R_out = p(1,end-3); % [Pa] - Pressure

    case {"subsonic_in_fixed_out"}

        % Fix outlet pressure, but let the rest of variables float:
        rho_R_in = 2*rho(1,end-2) - rho(1,end-3); % [kg/m3] - Density
        u_R_in = 2*u(1,end-2) - u(1,end-3); % [m/s] - Velocity
        p_R_in = 2*p_out - p(1,end-2); % [Pa] - Pressure

        % Get outer node variables (assuming transmissive conditions):
        rho_R_out = rho_R_in; % [kg/m3] - Density
        u_R_out = u_R_in; % [m/s] - Velocity
        p_R_out = p_R_in; % [Pa] - Pressure

    otherwise
        error("Boundary conditions are not valid! Try another configuration or define new ones.")
end

% Get outlet/right approximated gas properties:
Rg_R_in = Rg(1,end-1); % [J/(kgK)] - Gas constant (inner halo node)
Rg_R_out = Rg(1,end); % [J/(kgK)] - Gas constant (outer halo node)
gamma_R_in = gamma(1,end-1); % [] - Gas adiabatic coefficient (inner halo node)
gamma_R_out = gamma(1,end); % [] - Gas adiabatic coefficient (outer halo node)
cv_R_in = getcv(Rg_R_in,gamma_R_in); % [J/(kgK)] - Gas specific heat at constant volume (inner halo node)
cv_R_out = getcv(Rg_R_out,gamma_R_out); % [J/(kgK)] - Gas specific heat at constant volume (outer halo node)

% Compute corresponding temperatures T:
T_R_in = T_from_state_equation(rho_R_in,Rg_R_in,p_R_in); % [K] - Temperature (inner halo node)
T_R_out = T_from_state_equation(rho_R_out,Rg_R_out,p_R_out); % [K] - Temperature (outer halo node)

% Compute corresponding internal specific energies e:
e_R_in = cv_R_in*T_R_in; % [m2/s2] - Specific energy (inner halo node)
e_R_out = cv_R_out*T_R_out; % [m2/s2] - Specific energy (outer halo node)

end

% ----------------------------------------------------------------------- %