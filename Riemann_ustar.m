% ----------------------------------------------------------------------- %
% Function that gets the velocity at the Star Region u* of a Riemann
% problem, given left and right data and Star Region pressure p*.
% ----------------------------------------------------------------------- %
% Input:
% - gamma: gas ratio of specific heats.
% - WL: vector of left primitive variables [rhoL,uL,pL]^T [kg/m3,m/s,Pa].
% - WR: vector of right primitive variables [rhoR,uR,pR]^T [kg/m3,m/s,Pa].
% - pstar: pressure solution, i.e. pressure at the Star Region [Pa].
% ----------------------------------------------------------------------- %
% Output:
% - ustar: velocity solution, i.e. velocity at the Star Region [m/s].
% ----------------------------------------------------------------------- %

function ustar = Riemann_ustar(gamma,WL,WR,pstar)

% Getting left and right velocities:
uL = WL(2); % [m/s]
uR = WR(2); % [m/s]

% Getting left and right Riemann's algebraic equation results for p*:
[~,~,fL,fR] = Riemann_f_p(gamma,WL,WR,pstar);

% Computing u*:
ustar = 0.5*(uL + uR) + 0.5*(fR - fL); % [m/s]

end

% ----------------------------------------------------------------------- %