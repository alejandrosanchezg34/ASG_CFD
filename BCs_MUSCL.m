% ----------------------------------------------------------------------- %
% This function imposes the desired boundary conditions by adapting the
% flow-field variable vectors involved in the problem.
% This is expected to be implemented with the 1D or quasi-1D MUSCL-Hancock
% solver.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density vector [kg/m3] [1 x N+4].
% - u: velocity vector [m/s] [1 x N+4].
% - e: specific energy vector [m2/s2] [1 x N+4].
% - T: temperature vector [K] [1 x N+4].
% - p: pressure vector [Pa] [1 x N+4].
% - Rg: gas constant [J/(kgK)] [1 x N+4].
% - gamma: gas adiabatic coefficient [] [1 x N+4].
% - rho_in: if needed, inlet density to be fixed [kg/m3].
% - u_in: if needed, inlet velocity to be fixed [m/s].
% - p_in: if needed, inlet pressure to be fixed [Pa].
% - p_out: if needed, outlet pressure to be fixed [Pa].
% - BCs: string specifying the type of boundary conditions to impose.
% ----------------------------------------------------------------------- %
% Output:
% - rho: completed density vector [kg/m3] [1 x N+4].
% - u: completed velocity vector [m/s] [1 x N+4].
% - e: completed specific energy vector [m2/s2] [1 x N+4].
% - T: completed temperature vector [K] [1 x N+4].
% - p: completed pressure vector [Pa] [1 x N+4].
% ----------------------------------------------------------------------- %

function [rho,u,e,T,p] = BCs_MUSCL(rho,u,e,T,p,Rg,gamma,rho_in,u_in,p_in,p_out,BCs)

% Impose boundary conditions:

% --- Get inlet values
[rho_L_out,u_L_out,e_L_out,T_L_out,p_L_out, ...
    rho_L_in,u_L_in,e_L_in,T_L_in,p_L_in] = BCs_MUSCL_in(rho,u,p,Rg,gamma,rho_in,u_in,p_in,BCs);

% --- Get outlet values
[rho_R_in,u_R_in,e_R_in,T_R_in,p_R_in, ...
    rho_R_out,u_R_out,e_R_out,T_R_out,p_R_out] = BCs_MUSCL_out(rho,u,p,Rg,gamma,p_out,BCs);

% --- Complete the flow-field variable vectors
rho(1,1) = rho_L_out; rho(1,2) = rho_L_in; rho(1,end-1) = rho_R_in; rho(1,end) = rho_R_out; % [kg/m3] - Density
u(1,1) = u_L_out; u(1,2) = u_L_in; u(1,end-1) = u_R_in; u(1,end) = u_R_out; % [m/s] - Velocity
e(1,1) = e_L_out; e(1,2) = e_L_in; e(1,end-1) = e_R_in; e(1,end) = e_R_out; % [m2/s2] - Internal specific energy
T(1,1) = T_L_out; T(1,2) = T_L_in; T(1,end-1) = T_R_in; T(1,end) = T_R_out; % [K] - Temperature
p(1,1) = p_L_out; p(1,2) = p_L_in; p(1,end-1) = p_R_in; p(1,end) = p_R_out; % [Pa] - Pressure

end

% ----------------------------------------------------------------------- %