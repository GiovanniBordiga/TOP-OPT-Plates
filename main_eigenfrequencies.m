%% MAIN SCRIPT FOR EIGENFREQUENCY OPTIMIZATION
% This code implements an optimization procedure to find, given a certain
% amount of available material, the mass distribution of a plate which 
% maximizes its foundamental (lowest) eigenfrequency.
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
nelx = 100; nely = 50;      % number of plate elements 
dims.width = 1; dims.height = 1; dims.thickness = 1;            % element's dimensions
material.E = 2.1e+11; material.v = 0.3; material.rho = 7800;    % material properties
element = FE('MB4', dims, material);                            % build the finite element
FrVol = 0.3;                % volume fraction at the optimum condition
x = ones(nely, nelx)*FrVol; % set uniform intial density

%% PROBLEM SELECTION
problem = Problem(nelx, nely, element, 'a'); % list of problems in "FEM/Problem"
nModes = 5;                 % number of eigenmodes to compute

%% INITIALIZE NUMERICAL VARIABLES
PenK = 3;                   % stiffness penalization factor used in the SIMP model
PenM = 3;                   % mass penalization factor used in the SIMP model
RaFil = 2;                  % filter radius
move = 0.1;                 % limit to the change of 'x' (optimum)
SF = 0.1;                   % stabilization factor (optimum)

%% OPTIMIZATION CYCLE
tol = 1e-3;                 % tolerance for convergence criteria
change = 1;                 % density change in the plates (convergence)
changes = [];               % history of the density change (plot)
eigenFs = [];               % history of the eigenfrequencies (plot)
maxiter = 15;               % maximum number of iterations (convergence)
iter = 0;                   % iteration counter
while change > tol && iter < maxiter
    %% optimize
    [eigenF, eigenM] = eigenFM(problem, element,...
        x, PenK, PenM, nModes);                             % solve the eigenvalues problem
    dF = getFSensitivity(nelx, nely, element, x,...
        PenK, PenM, eigenF, eigenM);                        % sensitivity analysis
    dF = filterSensitivity(nelx, nely, x, dF, RaFil);           % apply sensitivity filter
    xnew = OC(nelx, nely, element, x, FrVol, -dF, move, SF);    % get new densities
    change = max(max(abs(xnew-x)));
    x = xnew;               % update densities
    iter = iter + 1;
    %% display results
    disp(['Iter: ' sprintf('%i', iter) ', Obj: ' sprintf('%.3f', min(eigenF))...
        ', Vol. frac.: ' sprintf('%.3f', sum(sum(x))/(nelx*nely))]);
    eigenFs = cat(2, eigenFs, min(eigenF));
    changes = cat(2, changes, change);
    plotConvergence(1:iter, eigenFs, 'f');
    plotConvergence(1:iter, changes, 'x');
    plotDesign(x);
end