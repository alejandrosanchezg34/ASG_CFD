% ----------------------------------------------------------------------- %
% Function that, given the flow-field variables, returns the corresponding
% source terms vector J.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - p: pressure [Pa].
% - dAdx: cross-section area derivative.
% * May be given as vectors [1 x N].
% ----------------------------------------------------------------------- %
% Output:
% - J: source terms vector ([3 x 1] or [3 x N]).
% ----------------------------------------------------------------------- %

function J = vars2J(p,dAdx)

% Number of nodes being evaluated:
N = size(p,2);

% Arrange J vector:
J = [zeros(1,N);
     p.*dAdx;
     zeros(1,N)];

end

% ----------------------------------------------------------------------- %