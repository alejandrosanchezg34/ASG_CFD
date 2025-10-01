% ----------------------------------------------------------------------- %
% Function that computes the slope vector used to reconstruct the data cell
% average values into piece-wise linear functions. This is done at a given
% cell 'i'.
% Function valid for 1D problems.
% Note that the function is able to implement different TVD methods, based
% on limited slopes or slope limiters.
% ----------------------------------------------------------------------- %
% Input:
% - i: ID of the cell being evaluated.
% - U: solution vectors matrix at the evaluated time instant [3 x N].
% - w: constant parameter (can go from -1 to 1).
% ----------------------------------------------------------------------- %
% Output:
% - di: slope vector at cell 'i' [3 x 1].
% ----------------------------------------------------------------------- %

function di = slope_vector(i,U,w)

% Select TVD method among the following:
% - "limited_slopes"
% - "slope_limiter_SUPERBEE"
% - "slope_limiter_MINBEE"
% - "none"
TVD = "none";

% Local left and right intercell slope vectors:
diL = U(:,i) - U(:,i-1); % Left intercell
diR = U(:,i+1) - U(:,i); % Right intercell

% Slope vector at cell 'i':
di = 0.5*(1+w)*diL + 0.5*(1-w)*diR;

% Apply TVD method:
switch TVD
    case "limited_slopes"
        
        % --- Limited Slopes

        % Select parameter beta:
        beta = 1; % MINBEE flux limiter
        % beta = 2; % SUPERBEE flux limiter
        
        % For each slope vector component:
        for j = 1:length(di)
            if diR(j) > 0
                cases = [0, min(beta*diL(j),diR(j)), min(diL(j),beta*diR(j))];
                di(j) = max(cases,[],"all");
        
            elseif diR(j) < 0
                cases = [0, max(beta*diL(j),diR(j)), max(diL(j),beta*diR(j))];
                di(j) = min(cases,[],"all");
        
            end
        end

    case "slope_limiter_SUPERBEE"
        
        % --- SUPERBEE Slope Limiter
        xi = zeros(length(di),1);

        % Define beta parameter at 'i+1/2':
        betaR = 1;

        % For each slope limiter component:
        for j = 1:length(di)
            r = diL(j)/diR(j); % Slopes ratio

            if (r <= 0)
                xi(j) = 0;

            elseif (r > 0 && r <= 0.5)
                xi(j) = 2*r;

            elseif (r > 0.5 && r <= 1)
                xi(j) = 1;

            else
                xiR = 2*betaR/(1-w+(1+w)*r);
                xi(j) = min([r,xiR,2]);

            end
        end

        % Correct slope vector:
        di = xi.*di;

    case "slope_limiter_MINBEE"
        
        % --- MINBEE Slope Limiter
        xi = zeros(length(di),1);

        % Define beta parameter at 'i+1/2':
        betaR = 1;

        % For each slope limiter component:
        for j = 1:length(di)
            r = diL(j)/diR(j); % Slopes ratio

            if (r <= 0)
                xi(j) = 0;

            elseif (r > 0 && r <= 1)
                xi(j) = r;

            else
                xiR = 2*betaR/(1-w+(1+w)*r);
                xi(j) = min([1,xiR]);

            end
        end

        % Correct slope vector:
        di = xi.*di;

    case "none"
        
        % Do not modify slope vector

    otherwise
        error("The TVD method specified in 'slope_vector.m' function is not available");

end

end

% ----------------------------------------------------------------------- %