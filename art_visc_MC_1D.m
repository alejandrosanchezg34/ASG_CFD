% ----------------------------------------------------------------------- %
% Function that computes MacCormack's artificial viscosity to be added to
% the Euler equations, according to: Anderson. "CFD: The Basics with
% Applications", 1995.
% Valid function for 1D and quasi-1D problems.
% ----------------------------------------------------------------------- %
% Input:
% - p: pressure vector [Pa] [1 x N].
% - U: solution vector [3 x N].
% - Cx: artificial viscosity coefficient.
% - i: index associated to the node being evaluated.
% ----------------------------------------------------------------------- %
% Output:
% - Si: MacCormack's artificial viscosity vector at given node 'i'.
% ----------------------------------------------------------------------- %

function Si = art_visc_MC_1D(p,U,Cx,i)

% Obtain artificial viscosity vector at given node 'i':
Si = Cx*abs(p(i+1) - 2*p(i) + p(i-1))/(p(i+1) + 2*p(i) + p(i-1))*(U(:,i+1) - 2*U(:,i) + U(:,i-1));

end

% ----------------------------------------------------------------------- %