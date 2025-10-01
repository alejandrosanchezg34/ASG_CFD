% ----------------------------------------------------------------------- %
% This function provides the exact solution of the Riemann problem at a
% given point in space 'x' and at a given time 't'. The solution is
% provided as a vector of primitive variables: density, velocity and
% pressure.
% ----------------------------------------------------------------------- %
% Input:
% - x: position in space where to evaluate the problem, between -L and L
% [m].
% - t: time instant when to evaluate the problem, starting from 0 [s].
% - gamma: gas ratio of specific heats.
% - WL: vector of left primitive variables [rhoL,uL,pL]^T [kg/m3,m/s,Pa].
% - WR: vector of right primitive variables [rhoR,uR,pR]^T [kg/m3,m/s,Pa].
% ----------------------------------------------------------------------- %
% Output:
% - WSol: solution vector of primitive variables [rho,u,p]^T [kg/m3,m/s,Pa].
% ----------------------------------------------------------------------- %

function WSol = Riemann_solution(x,t,gamma,WL,WR)

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
% ----------------------------------------------------------------------- %

% Obtain p* and u*:
pstar = Riemann_pstar(gamma,WL,WR); % [Pa]
ustar = Riemann_ustar(gamma,WL,WR,pstar); % [m/s]

% Compute speed (term needed to sample the solution):
S = x/t; % [m/s]

% Two possible cases:

if (S <= ustar)
    % Sampling at left side of contact:
    
    if (pstar > pL)
        % Left shock wave:
        SL = uL - aL*((gamma+1)/(2*gamma)*pstar/pL + (gamma-1)/(2*gamma))^0.5; % [m/s] - Shock speed
        
        if (S <= SL) % Before shock            
            WSol = WL;

        else % After shock            
            rhoLstar = rhoL*((pstar/pL + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*pstar/pL + 1)); % [kg/m3] - Density
            WSol = [rhoLstar; ustar; pstar];       

        end

    else
        % Left rarefaction wave:
        rhoLstar = rhoL*(pstar/pL)^(1/gamma); % [kg/m3] - Density
        aLstar = sqrt(gamma*pstar/rhoLstar); % [m/s] - Sound speed

        SHL = uL - aL; % [m/s] - Head speed of rarefaction wave
        STL = ustar - aLstar; % [m/s] - Tail speed of rarefaction wave

        if (S <= SHL) % Before rarefaction
            WSol = WL;

        elseif (S > SHL && S < STL) % Inside rarefaction
            rho_fan = rhoL*(2/(gamma+1) + (gamma-1)/((gamma+1)*aL)*(uL - S))^(2/(gamma-1)); % [kg/m3] - Density
            u_fan = 2/(gamma+1)*(aL + (gamma-1)/2*uL + S); % [m/s] - Velocity
            p_fan = pL*(2/(gamma+1) + (gamma-1)/((gamma+1)*aL)*(uL - S))^(2*gamma/(gamma-1)); % [Pa] - Pressure
            WSol = [rho_fan; u_fan; p_fan];

        else % After rarefaction
            WSol = [rhoLstar; ustar; pstar];

        end
    end
else
    % Sampling at right side of contact:

    if (pstar > pR)
        % Right shock wave:
        SR = uR + aR*((gamma+1)/(2*gamma)*pstar/pR + (gamma-1)/(2*gamma))^0.5; % [m/s] - Shock speed
        
        if (S >= SR) % Before shock            
            WSol = WR;

        else % After shock
            rhoRstar = rhoR*((pstar/pR + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*pstar/pR + 1)); % [kg/m3] - Density
            WSol = [rhoRstar; ustar; pstar];       

        end

    else
        % Right rarefaction wave:
        rhoRstar = rhoR*(pstar/pR)^(1/gamma); % [kg/m3] - Density
        aRstar = sqrt(gamma*pstar/rhoRstar); % [m/s] - Sound speed

        SHR = uR + aR; % [m/s] - Head speed of rarefaction wave
        STR = ustar + aRstar; % [m/s] - Tail speed of rarefaction wave

        if (S >= SHR) % Before rarefaction
            WSol = WR;

        elseif (S < SHR && S > STR) % Inside rarefaction
            rho_fan = rhoR*(2/(gamma+1) - (gamma-1)/((gamma+1)*aR)*(uR - S))^(2/(gamma-1)); % [kg/m3] - Density
            u_fan = 2/(gamma+1)*(-aR + (gamma-1)/2*uR + S); % [m/s] - Velocity
            p_fan = pR*(2/(gamma+1) - (gamma-1)/((gamma+1)*aR)*(uR - S))^(2*gamma/(gamma-1)); % [Pa] - Pressure
            WSol = [rho_fan; u_fan; p_fan];

        else % After rarefaction
            WSol = [rhoRstar; ustar; pstar];

        end        
    end
end

end

% ----------------------------------------------------------------------- %