% ----------------------------------------------------------------------- %
% This function implements the HLLC fluxes approximate Riemann solver for a
% 1D Riemann problem. Given the sets of left and right data, the resulting 
% HLLC intercell flux F (at i+1/2) is returned.
% ----------------------------------------------------------------------- %
% Input:
% - WL: vector of left primitive variables [rhoL,uL,pL]^T [kg/m3,m/s,Pa].
% - WR: vector of right primitive variables [rhoR,uR,pR]^T [kg/m3,m/s,Pa].
% - Rg: gas constant [J/(kgK)].
% - gamma: gas ratio of specific heats.
% - opt: output options. If "SK" is written, only the non-linear wave
% speeds are returned. If nothing is written, the whole solution is
% computed.
% ----------------------------------------------------------------------- %
% Output:
% - Fint: HLLC intercell flux vector (at i+1/2) [3 x 1].
% - SL: left non-linear wave speed [m/s].
% - SR: right non-linear wave speed [m/s].
% - Uint: HLLC intercell solution vector (at i+1/2) [3 x 1].
% ----------------------------------------------------------------------- %

function [Fint,SL,SR,Uint] = HLLC_solver(WL,WR,Rg,gamma,opt)

if nargin == 4
    opt = "all";
end

% ----------------------------------------------------------------------- %

% Compute gas specific heat at constant volume:
cv = Rg/(gamma - 1); % [J/(kgK)]

% Obtaining left variables:
rhoL = WL(1); % [kg/m3]
uL = WL(2); % [m/s]
pL = WL(3); % [Pa]
TL = pL/(rhoL*Rg); % [K]
eL = cv*TL; % [m2/s2]

% Obtaining right variables:
rhoR = WR(1); % [kg/m3]
uR = WR(2); % [m/s]
pR = WR(3); % [Pa]
TR = pR/(rhoR*Rg); % [K]
eR = cv*TR; % [m2/s2]

% Sound speeds:
aL = sqrt(gamma*pL/rhoL); % [m/s]
aR = sqrt(gamma*pR/rhoR); % [m/s]

% ----------------------------------------------------------------------- %

% Density and sound speed averages:
rho_avg = 0.5*(rhoL + rhoR); % [kg/m3]
a_avg = 0.5*(aL + aR); % [m/s]

% Star Region pressure estimate (based on PVRS):
p_pvrs = 0.5*(pL + pR) - 0.5*(uR - uL)*rho_avg*a_avg; % [Pa]
pstar = max(0,p_pvrs); % [Pa]

% ----------------------------------------------------------------------- %

% Wave speed estimates:
qL = qK_fun(pL,pstar,gamma); % Left qK term
qR = qK_fun(pR,pstar,gamma); % Right qK term

SL = uL - aL*qL; % [m/s] - Left non-linear wave speed
SR = uR + aR*qR; % [m/s] - Right non-linear wave speed

if opt == "SK"
    Fint = NaN;
    Uint = NaN;
    return;
end

% Intermediate speed:
Sstar = (pR - pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR))/(rhoL*(SL - uL) - rhoR*(SR - uR)); % [m/s]

% ----------------------------------------------------------------------- %

% Left and right solution U and flux F vectors:
UL = vars2U(rhoL,uL,eL,1);
UR = vars2U(rhoR,uR,eR,1);
FL = vars2F(rhoL,uL,eL,pL,1);
FR = vars2F(rhoR,uR,eR,pR,1);

% HLLC intercell flux Fint and solution vector Uint (at local x/t = 0 line):
if (0 <= SL)
    Uint = UL;
    Fint = FL;

elseif (SL <= 0 && 0 <= Sstar)
    [Uint,Fint] = U_F_starK_fun(UL,FL,SL,Sstar,Rg,gamma);

elseif (Sstar <= 0 && 0 <= SR)
    [Uint,Fint] = U_F_starK_fun(UR,FR,SR,Sstar,Rg,gamma);

elseif (0 >= SR)
    Uint = UR;
    Fint = FR;

end

end

% ----------------------------------------------------------------------- %

%% Other functions

% qK term used to estimate the non-linear wave speeds:
function qK = qK_fun(pK,pstar,gamma)

if pstar <= pK
    qK = 1;
else
    qK = (1 + ((gamma+1)/(2*gamma))*(pstar/pK - 1))^0.5;
end

end

% Star Region HLLC flux and solution vectors:
function [UstarK,FstarK] = U_F_starK_fun(UK,FK,SK,Sstar,Rg,gamma)

% Flow-field K variables:
[rhoK,uK,eK,~,pK] = U2vars(UK,1,Rg,gamma);

% Star region solution vector UstarK:
UstarK = rhoK*(SK - uK)/(SK - Sstar)*[1;
                                      Sstar;
                                      (eK + uK^2/2) + (Sstar - uK)*(Sstar + pK/(rhoK*(SK - uK)))];

% Star region flux vector FstarK:
FstarK = FK + SK*(UstarK - UK);

end

% ----------------------------------------------------------------------- %