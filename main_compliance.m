import FEM.*
import opt.*

%% INITIALIZE GEOMETRY & MATERIAL
nelx = 40; nely = 10;               % number of plate elements 
dims.width = 1; dims.height = 1; dims.thickness = 1; % element's dimensions
element = 'ACM';                    % finite element type
material.E = 1; material.v = 0.3;   % material properties
FrVol = 0.3;                % volume fraction at the optimum condition

%% INITIALIZE DESIGN VARIABLE
CoPen = 3;                  % penalization coefficient used in the SIMP model
x = ones(nely, nelx)*FrVol; % intial value of the density field

%% OPTIMIZATION CYCLE
change = 1;                 % maximum density change in the plates (convergence)
maxiter = 5;                % maximum number of iteration (convergence)
iter = 0;                   % iteration counter
Ke = getK(element, dims, material);
while change > 1e-3 && iter < maxiter
    % TODO optimize!
    U = FEM(nelx, nely, element, dims, material, x, CoPen);
    [dC, C] = getSensitivity(nelx, nely, x, CoPen, Ke, U);
    xnew = OC(nelx, nely, x, FrVol, dC);
    change = max(max(abs(xnew-x)));
    x = xnew;
    iter = iter + 1;
    disp(['Iter: ' sprintf('%i', iter) ', Obj: ' sprintf('%.3f', C)...
        ', Vol. frac.: ' sprintf('%.3f', sum(sum(x))/(nelx*nely))]);
end