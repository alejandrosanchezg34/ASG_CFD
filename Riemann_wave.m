% ----------------------------------------------------------------------- %
% This function provides the non-linear wave speeds present at the solution
% of a Riemann problem of both left and right waves. In the case of
% rarefaction waves, the head speed is returned.
% ----------------------------------------------------------------------- %
% Input:
% - gamma: gas ratio of specific heats.
% - WL: vector of left primitive variables [rhoL,uL,pL]^T [kg/m3,m/s,Pa].
% - WR: vector of right primitive variables [rhoR,uR,pR]^T [kg/m3,m/s,Pa].
% ----------------------------------------------------------------------- %
% Output:
% - SL: wave speed of the left non-linear wave [m/s].
% - SR: wave speed of the right non-linear wave [m/s].
% ----------------------------------------------------------------------- %

function [SL,SR] = Riemann_wave(gamma,WL,WR)

% ----------------------------------------------------------------------- %
% Obtaining left variables:
rhoL = WL(1); % [kg/m3]
uL = WL(2); % [m/s]
pL = WL(3); % [Pa]

% Obtaining right variables:
rhoR = WR(1); % [kg/m3]
uR = WR(2); % [m/s]
pR = WR(3); % [Pa]

% Sound speeds:
aL = sqrt(gamma*pL/rhoL); % [m/s]
aR = sqrt(gamma*pR/rhoR); % [m/s]

% Obtain p*:
pstar = Riemann_pstar(gamma,WL,WR); % [Pa]
% ----------------------------------------------------------------------- %

% Sampling at left side of contact:
if (pstar > pL)    
    % Left shock wave:
    SL = uL - aL*((gamma+1)/(2*gamma)*pstar/pL + (gamma-1)/(2*gamma))^0.5; % [m/s] - Shock speed
else
    % Left rarefaction wave:
    SL = uL - aL; % [m/s] - Head speed of rarefaction wave
end

% Sampling at right side of contact:
if (pstar > pR)
    % Right shock wave:
    SR = uR + aR*((gamma+1)/(2*gamma)*pstar/pR + (gamma-1)/(2*gamma))^0.5; % [m/s] - Shock speed
else
    % Right rarefaction wave:
    SR = uR + aR; % [m/s] - Head speed of rarefaction wave
end

end

% ----------------------------------------------------------------------- %