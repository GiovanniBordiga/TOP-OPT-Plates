import FEM.*
import opt.*

%% PROBLEM SELECTION
problem = 'a';

%% INITIALIZE GEOMETRY, MATERIAL, DESIGN VARIABLE
nelx = 40; nely = 10;               % number of plate elements 
dims.width = 1; dims.height = 1; dims.thickness = 1; % element's dimensions
element = 'ACM';                    % finite element type
material.E = 1; material.v = 0.3;   % material properties
FrVol = 0.3;                % volume fraction at the optimum condition
x = ones(nely, nelx)*FrVol; % set uniform intial density

%% INITIALIZE NUMERICAL VARIABLES
CoPen = 3;                  % penalization coefficient used in the SIMP model
RaFil = 2;                  % filter radius

%% OPTIMIZATION CYCLE
change = 1;                 % maximum density change in the plates (convergence)
maxiter = 5;                % maximum number of iteration (convergence)
iter = 0;                   % iteration counter
[Kf, Ks] = getK(element, dims, material); Ke = Kf + Ks;
while change > 1e-3 && iter < maxiter
    U = FEM(problem, nelx, nely, element, dims, material, x, CoPen); % solve FEM
    [dC, C] = getSensitivity(nelx, nely, x, CoPen, Ke, U);  % sensitivity analysis
    dC = filterSensitivity(nelx, nely, x, dC, RaFil);       % apply sensitivity filter
    xnew = OC(nelx, nely, x, FrVol, dC);                    % get new densities
    change = max(max(abs(xnew-x)));
    x = xnew;           % update densities
    iter = iter + 1;
    disp(['Iter: ' sprintf('%i', iter) ', Obj: ' sprintf('%.3f', C)...
        ', Vol. frac.: ' sprintf('%.3f', sum(sum(x))/(nelx*nely))]);
end