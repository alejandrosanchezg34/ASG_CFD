% ----------------------------------------------------------------------- %
% This function imposes the desired boundary conditions by adapting the
% flow-field variable vectors involved in the problem.
% This is expected to be implemented with the 1D or quasi-1D Godunov
% solver.
% ----------------------------------------------------------------------- %
% Input:
% - rho: density vector [kg/m3] [1 x N+2].
% - u: velocity vector [m/s] [1 x N+2].
% - e: specific energy vector [m2/s2] [1 x N+2].
% - T: temperature vector [K] [1 x N+2].
% - p: pressure vector [Pa] [1 x N+2].
% - Rg: gas constant [J/(kgK)] [1 x N+2].
% - gamma: gas adiabatic coefficient [] [1 x N+2].
% - rho_in: if needed, inlet density to be fixed [kg/m3].
% - u_in: if needed, inlet velocity to be fixed [m/s].
% - p_in: if needed, inlet pressure to be fixed [Pa].
% - p_out: if needed, outlet pressure to be fixed [Pa].
% - BCs: string specifying the type of boundary conditions to impose.
% - Riemann_solver: Riemann problem solver to be used ("HLLC" or "exact").
% ----------------------------------------------------------------------- %
% Output:
% - rho: completed density vector [kg/m3] [1 x N+2].
% - u: completed velocity vector [m/s] [1 x N+2].
% - e: completed specific energy vector [m2/s2] [1 x N+2].
% - T: completed temperature vector [K] [1 x N+2].
% - p: completed pressure vector [Pa] [1 x N+2].
% ----------------------------------------------------------------------- %

function [rho,u,e,T,p] = BCs_Godunov(rho,u,e,T,p,Rg,gamma,rho_in,u_in,p_in,p_out,BCs,Riemann_solver)

% Impose boundary conditions:

% --- Get inlet values
[rho_L,u_L,e_L,T_L,p_L] = BCs_Godunov_in(rho,u,p,Rg,gamma,rho_in,u_in,p_in,BCs,Riemann_solver);

% --- Get outlet values
[rho_R,u_R,e_R,T_R,p_R] = BCs_Godunov_out(rho,u,p,Rg,gamma,p_out,BCs,Riemann_solver);

% --- Complete the flow-field variable vectors
rho(1,1) = rho_L; rho(1,end) = rho_R; % [kg/m3] - Density
u(1,1) = u_L; u(1,end) = u_R; % [m/s] - Velocity
e(1,1) = e_L; e(1,end) = e_R; % [m2/s2] - Internal specific energy
T(1,1) = T_L; T(1,end) = T_R; % [K] - Temperature
p(1,1) = p_L; p(1,end) = p_R; % [Pa] - Pressure

end

% ----------------------------------------------------------------------- %