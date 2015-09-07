function xnew = PG(nelx, nely, x, FrVol, dF)
%PG Summary of this function goes here
%   Detailed explanation goes here

stepsize = 1;
l = sum(sum(dF))/(nelx*nely);
% stepsize = (FrVol*nelx*nely - sum(sum(x)))/(sum(sum(-dF + l)));
xnew = max(0.001, min(1, x + stepsize*(-dF + l)));
end

