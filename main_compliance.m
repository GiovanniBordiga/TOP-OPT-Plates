%% MAIN SCRIPT FOR COMPLIANCE OPTIMIZATION
% This code implements an optimization procedure to find, given a certain
% amount of available material, the mass distribution of a plate which 
% minimizes its compliance (work of loads).
% Dimensions and material of the plate can be specified in the following
% initialization section. Only rectangular plates are considered.
% Three types of finite element are supported: two Kirchhoff elements
% (ACM and BMF) and one Mindlin element (here called MB4). To use the
% desired finite element just plug the corresponding string in the
% initialization of the 'FE' object.
% Boundary conditions (constraints and loads) are specified in the
% 'Problem' class and initialized through its object 'problem' passing the
% appropriate string ('problemId') to the constructor. For a full list of
% available problems, see the 'Problem' class in the 'FEM' package.
% To set the volume constraint of the optimization problem, use the volume
% fraction variable 'FrVol'.
%%

import FEM.*
import opt.*
import plot.*

%% INITIALIZE GEOMETRY, MATERIAL, DESIGN VARIABLE
nelx = 100; nely = 100;      % number of plate elements along the two axes
dims.width = 1; dims.height = 1; dims.thickness = 1;    % element's dimensions
material.E = 1000; material.v = 0.3; material.rho = 1;  % material properties
element = FE('MB4', dims, material);                    % build the finite element
FrVol = 0.3;                % volume fraction at the optimum condition
x = ones(nely, nelx)*FrVol; % set uniform intial density

%% PROBLEM SELECTION
problem = Problem(nelx, nely, element, 'e'); % list of problems in "FEM/Problem"

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
while change > tol && iter < maxiter
    %% optimize
    U = FEM(problem, element, x, CoPen);                            % solve FEM
    [dC, C] = getCSensitivity(nelx, nely, element, x, CoPen, U);    % sensitivity analysis
    dC = filterSensitivity(nelx, nely, x, dC, RaFil);               % apply sensitivity filter
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
plotDeformed(nelx, nely, element, x.^CoPen, U);
