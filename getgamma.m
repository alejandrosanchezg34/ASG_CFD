% ----------------------------------------------------------------------- %
% This function, for given values of T and p, returns the corresponding gas
% adiabatic coefficient 'gamma'. Three different gas models are implemented:
% - Constant: assumes perfect gases, i.e. constant properties and
% composition with T and p.
% - Frozen flow: under the assumption of ideal gases, assumes the gas
% properties to vary with T, but the composition to remain constant.
% - Shifting flow: under the assumption of ideal gases, assumes both the
% gas properties and the composition to vary with T.
% The last two models are implemented with HGS.
% ----------------------------------------------------------------------- %
% Input:
% - gas_model: string variable that specifies the gas model selected. It
% can be either "constant", "frozen_flow" or "shifting_flow".
% - gamma_ct: gas adiabatic coefficient to be used by the constant gas
% model [].
% - species: cell array with the names (strings) of the species involved in
% the gas. It is not needed for the "constant" gas model.
% - n_0: initial/chamber gas composition [mols] (given as a [species x 1]
% vector).
% - T: temperature [K] (may be given as a [1 x N] vector).
% - p: pressure [Pa] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %
% Output:
% - gamma: gas adiabatic coefficient [] (may be given as a [1 x N] vector).
% ----------------------------------------------------------------------- %

function gamma = getgamma(gas_model,gamma_ct,species,n_0,T,p)

if (nargin == 2 && gas_model == "constant")
    species = NaN; % Not needed
    n_0 = NaN; % Not needed
    T = NaN; % Not needed
    p = NaN; % Not needed
end

% Number of grid points:
N = size(T,2);

% Compute gamma according to the gas model:
gamma = zeros(1,N); % Initialize vector

switch gas_model
    case "constant"
        % Set the constant value for gamma at each node (no dependence on T):
        gamma(1,:) = gamma_ct; % []

    case "frozen_flow"
        % Call HGS for every temperature and pressure values:
        for i = 1:N
            % Avoid consulting HGS if 'T = NaN' and/or 'p = NaN':
            if (isnan(T(1,i)) || isnan(p(1,i)))
                gamma(1,i) = NaN;
                continue;
            end

            % Get 'gamma':
            gamma(1,i) = HGSprop(species,n_0,T(1,i),p(1,i)*1E-5,'gamma'); % []
        end

    case "shifting_flow"
        % Call HGS for every temperature and pressure values:
        for i = 1:N
            % Avoid consulting HGS if 'T = NaN' and/or 'p = NaN':
            if (isnan(T(1,i)) || isnan(p(1,i)))
                gamma(1,i) = NaN;
                continue;
            end

            % Calculate the equilibrium state of the mixture at given
            % temperature and pressure:
            [~,n_eq] = HGSeq(species,n_0,T(1,i),p(1,i)*1E-5); % [mols] - Equilibrium composition
            
            % Get 'gamma' for the equilibrium composition:
            gamma(1,i) = HGSprop(species,n_eq,T(1,i),p(1,i)*1E-5,'gamma'); % []
        end

    otherwise
        error("The gas model selected to compute gamma is wrong or not available!");
end

end

% ----------------------------------------------------------------------- %