import FEM.*
import opt.*
import plot.*

%% INITIALIZE GEOMETRY, MATERIAL, DESIGN VARIABLE
nelx = 100; nely = 50;       % number of plate elements 
dims.width = 1; dims.height = 1; dims.thickness = 1; % element's dimensions
element = FE('MB4', dims);  % build the finite element
material.E = 2.1e+11; material.v = 0.3; material.rho = 7800;   % material properties
FrVol = 0.3;                % volume fraction at the optimum condition
x = ones(nely, nelx)*FrVol; % set uniform intial density

%% PROBLEM SELECTION
problem = Problem(nelx, nely, element, 'a'); % list of problems in "FEM/Problem"
nModes = 1;                 % number of eigenmodes to compute

%% INITIALIZE NUMERICAL VARIABLES
PenK = 3;                   % stiffness penalization factor used in the SIMP model
PenM = 3;                   % mass penalization factor used in the SIMP model
RaFil = 2;                  % filter radius

%% OPTIMIZATION CYCLE
tol = 1e-3;                 % tolerance for convergence criteria
change = 1;                 % density change in the plates (convergence)
changes = [];               % history of the density change (plot)
eigenFs = [];               % history of the eigenfrequencies (plot)
maxiter = 5;               % maximum number of iterations (convergence)
iter = 0;                   % iteration counter
Ke = getK(element, material);
Me = getM(element, material);
while change > tol && iter < maxiter
    %% optimize
    [eigenF, eigenM] = eigenFM(problem, nelx, nely, element, material,...
        x, PenK, PenM, nModes);                             % solve the eigenvalues problem
    dF = getFSensitivity(nelx, nely, element, x,...
        PenK, PenM, Ke, Me, eigenF, eigenM);                % sensitivity analysis
    dF = filterSensitivity(nelx, nely, x, dF, RaFil);       % apply sensitivity filter
    xnew = OC(nelx, nely, x, FrVol, -dF);                   % get new densities
    change = max(max(abs(xnew-x)));
    x = xnew;               % update densities
    iter = iter + 1;
    %% display results
    disp(['Iter: ' sprintf('%i', iter) ', Obj: ' sprintf('%.3f', eigenF(end))...
        ', Vol. frac.: ' sprintf('%.3f', sum(sum(x))/(nelx*nely))]);
    eigenFs = cat(2, eigenFs, eigenF(end));
    changes = cat(2, changes, change);
    plotConvergence(1:iter, eigenFs, 'f');
    plotConvergence(1:iter, changes, 'x');
    plotDesign(x);
end