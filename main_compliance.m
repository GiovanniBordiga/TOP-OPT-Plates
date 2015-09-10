import FEM.*
import opt.*
import plot.*

%% INITIALIZE GEOMETRY, MATERIAL, DESIGN VARIABLE
nelx = 100; nely = 50;       % number of plate elements 
dims.width = 1; dims.height = 1; dims.thickness = 1; % element's dimensions
element = FE('ACM', dims);  % build the finite element
material.E = 1000; material.v = 0.3;   % material properties
FrVol = 0.3;                % volume fraction at the optimum condition
x = ones(nely, nelx)*FrVol; % set uniform intial density

%% PROBLEM SELECTION
problem = Problem(nelx, nely, element, 'a'); % list of problems in "FEM/Problem"

%% INITIALIZE NUMERICAL VARIABLES
CoPen = 3;                  % penalization coefficient used in the SIMP model
RaFil = 2;                  % filter radius
move = 0.2;                 % limit to the change of 'x' (optimum)
stepsize = 1;               % step factor along the gradient (optimum)

%% OPTIMIZATION CYCLE
tol = 1e-3;                 % tolerance for convergence criteria
change = 1;                 % density change in the plates (convergence)
changes = [];               % history of the density change (plot)
Cs = [];                    % history of the compliance (plot)
maxiter = 15;                % maximum number of iterations (convergence)
iter = 0;                   % iteration counter
Ke = getK(element, material);
while change > tol && iter < maxiter
    %% optimize
    U = FEM(problem, nelx, nely, element, material, x, CoPen); % solve FEM
    [dC, C] = getCSensitivity(nelx, nely, element, x, CoPen, Ke, U);  % sensitivity analysis
    dC = filterSensitivity(nelx, nely, x, dC, RaFil);       % apply sensitivity filter
    xnew = PG(nelx, nely, x, FrVol, dC, move, stepsize);    % get new densities
    change = max(max(abs(xnew-x)));
    x = xnew;               % update densities
    iter = iter + 1;
    %% display results
    disp(['Iter: ' sprintf('%i', iter) ', Obj: ' sprintf('%.3f', C)...
        ', Vol. frac.: ' sprintf('%.3f', sum(sum(x))/(nelx*nely))]);
    Cs = cat(2, Cs, C);
    changes = cat(2, changes, change);
    plotConvergence(1:iter, Cs, 'c');
    plotConvergence(1:iter, changes, 'x');
    plotDesign(x);
end

%% DISPLAY DEFORMED CONFIGURATION
plotDeformed(nelx, nely, element, x.^CoPen, Ke, U);
