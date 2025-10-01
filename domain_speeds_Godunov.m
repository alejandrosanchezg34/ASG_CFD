% ----------------------------------------------------------------------- %
% Function that computes all non-linear wave speeds present throughout the
% given domain. These wave speeds are saved and returned for further
% calculation of the required and stable time steps 'dt' in Godunov-type
% solver schemes.
% All vectors of primitive variables W throughout the domain are given as
% well.
% ----------------------------------------------------------------------- %
% Input:
% - U: solution vector [3 x Nt].
% - A: cross-section area distribution vector [1 x Nt].
% - Rg: gas constant [J/(kgK)] [1 x Nt].
% - gamma: gas adiabatic coefficient [1 x Nt].
% - Riemann_solver: Riemann problem solver to be used ("HLLC" or "exact").
% ----------------------------------------------------------------------- %
% Output:
% - SL: vector of speeds associated to the left non-linear waves [m/s]
% (from a 'i+1/2' Riemann problem) [Nt-1 x 1].
% - SR: vector of speeds associated to the right non-linear waves [m/s]
% (from a 'i+1/2' Riemann problem) [Nt-1 x 1].
% - W: vectors of primitive variables throughout the domain [3 x Nt].
% ----------------------------------------------------------------------- %

function [SL,SR,W] = domain_speeds_Godunov(U,A,Rg,gamma,Riemann_solver)

% Number of total nodes being evaluated (including additional ghost nodes):
Nt = size(U,2);

% Initialize wave speed and Riemann problem data vectors:
SL = zeros(Nt-1,1);
SR = zeros(Nt-1,1);
W = zeros(3,Nt);

% Find Riemann problem data throughout domain:
for i = 1:Nt
    W(:,i) = U2W(U(:,i),A(i),Rg(1,i),gamma(1,i));
end

% Save wave speeds found at 'i+1/2' Riemann problems:
for i = 1:(Nt-1)
    Rg_int = harmmean([Rg(1,i),Rg(1,i+1)]); % [J/(kgK)] - Intercell 'i+1/2' gas constant
    gamma_int = harmmean([gamma(1,i),gamma(1,i+1)]); % [] - Intercell 'i+1/2' gas adiabatic coefficient

    switch Riemann_solver
        case "HLLC"
            [~,SL(i),SR(i)] = HLLC_solver(W(:,i),W(:,i+1),Rg_int,gamma_int,"SK");

        case "exact"
            [SL(i),SR(i)] = Riemann_wave(gamma_int,W(:,i),W(:,i+1));
    end
end

end

% ----------------------------------------------------------------------- %