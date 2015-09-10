function xnew = PG(nelx, nely, x, FrVol, dF, move, stepsize)
% Implement a Projected Gradient method.
% 'nelx' and 'nely' are the number of elements along the two dimensions.
% 'x' is a nely-by-nelx matrix representing the density field on the plate.
% 'FrVol' is the desired volume fraction at the optimum condition.
% 'dF' is a nely-by-nelx matrix returned by the sensitivity analysis.

% WARNING: this code is written assuming unit element's area.
%          to minimize compliance pass dF.
%          to maximize eigenfrequency pass -dF.

l = sum(sum(dF))/(nelx*nely);
% stepsize = (FrVol*nelx*nely - sum(sum(x)))/(sum(sum(-dF + l)));
xnew = max(0.001, max(x-move, min(1, min(x+move, x + stepsize*(-dF + l)))));
end

