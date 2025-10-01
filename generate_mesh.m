% ----------------------------------------------------------------------- %
% Function that discretises the computational domain to be considered in
% the problems, by generating a mesh for a given number of nodes 'N'.
% The user can choose between a finite-difference mesh (to be used for the
% MacCormack technique) and a finite-volume mesh (to be used for the
% Godunov and MUSCL-Hancock schemes).
% ----------------------------------------------------------------------- %
% Input:
% - L: domain length [m].
% - N: number of grid points or number of cells.
% - method: either "MacCormack" or "Godunov".
% - x0: initial/left coordinate [m]. If it is not given, 'x0 = 0' is
% assumed.
% ----------------------------------------------------------------------- %
% Output:
% - x: vector of grid coordinates [m].
% - dx: cell width [m].
% ----------------------------------------------------------------------- %

function [x,dx] = generate_mesh(L,N,method,x0)

% If 'x0' is not given, it is assumed to be 'x0 = 0':
if nargin == 3
    x0 = 0; % [m]
end

% Vector of grid coordinates:
x = zeros(1,N);

% Select method:
switch method
    case "MacCormack"
        
        dx = L/(N-1); % [m] - Cell width

        % Generating mesh:
        for i = 1:N
            x(i) = (i-1)*dx + x0; % [m] - Grid coordinate
        end
        
    case "Godunov"
        
        dx = L/N; % [m] - Cell width

        % Generating mesh:
        for i = 1:N
            x(i) = (i-0.5)*dx + x0; % [m] - Grid coordinate
        end

    otherwise
        error("The method selected to generate the mesh is wrong!")
end

end

% ----------------------------------------------------------------------- %