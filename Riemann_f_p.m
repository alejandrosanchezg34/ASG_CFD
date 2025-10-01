% ----------------------------------------------------------------------- %
% Function that defines the algebraic equation as a function of pressure
% that is used to find the pressure in the Star Region p*, i.e. the
% solution for pressure of the Riemann problem. Its first derivative is
% also returned, as well as the left and right functions results.
% ----------------------------------------------------------------------- %
% Input:
% - gamma: gas ratio of specific heats.
% - WL: vector of left primitive variables [rhoL,uL,pL]^T [kg/m3,m/s,Pa].
% - WR: vector of right primitive variables [rhoR,uR,pR]^T [kg/m3,m/s,Pa].
% - p: pressure variable [Pa].
% ----------------------------------------------------------------------- %
% Output:
% - f: value associated to the algebraic equation for the given pressure.
% - fp: value associated to the first derivative of the algebraic equation
% for the given pressure.
% - fL: value associated to the left algebraic equation for the given
% pressure.
% - fR: value associated to the right algebraic equation for the given
% pressure.
% ----------------------------------------------------------------------- %

function [f,fp,fL,fR] = Riemann_f_p(gamma,WL,WR,p)

% ----------------------------------------------------------------------- %
% Obtaining left variables:
rhoL = WL(1); % [kg/m3]
uL = WL(2); % [m/s]
pL = WL(3); % [Pa]

% Obtaining right variables:
rhoR = WR(1); % [kg/m3]
uR = WR(2); % [m/s]
pR = WR(3); % [Pa]

% Velocity jump:
du = uR - uL; % [m/s]

% Sound speeds:
aL = sqrt(gamma*pL/rhoL); % [m/s]
aR = sqrt(gamma*pR/rhoR); % [m/s]

% Constants:
AL = 2/((gamma+1)*rhoL);
BL = ((gamma-1)/(gamma+1))*pL;
AR = 2/((gamma+1)*rhoR);
BR = ((gamma-1)/(gamma+1))*pR;
% ----------------------------------------------------------------------- %

% Left function:
if p > pL
    fL = (p - pL)*(AL/(p+BL))^0.5; % Left shock
else
    fL = 2*aL/(gamma-1)*((p/pL)^((gamma-1)/(2*gamma)) - 1); % Left rarefaction
end

% Right function:
if p > pR
    fR = (p - pR)*(AR/(p+BR))^0.5; % Right shock
else
    fR = 2*aR/(gamma-1)*((p/pR)^((gamma-1)/(2*gamma)) - 1); % Right rarefaction
end

% Algebraic function:
f = fL + fR + du;

% ----------------------------------------------------------------------- %

% Left function (1st derivative):
if p > pL
    fL_p = ((AL/(BL+p))^0.5)*(1 - (p - pL)/(2*(BL+p))); % Left shock
else
    fL_p = 1/(rhoL*aL)*(p/pL)^(-(gamma+1)/(2*gamma)); % Left rarefaction
end

% Right function (1st derivative):
if p > pR
    fR_p = ((AR/(BR+p))^0.5)*(1 - (p - pR)/(2*(BR+p))); % Right shock
else
    fR_p = 1/(rhoR*aR)*(p/pR)^(-(gamma+1)/(2*gamma)); % Right rarefaction
end

% First derivative of algebraic function:
fp = fL_p + fR_p;

end

% ----------------------------------------------------------------------- %