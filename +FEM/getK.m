function Ke = getK(element, material)
% Return the stiffness matrix of the given element.
% 'element' is a FE object.
% 'material' is a struct containing three attributes: E (Young's modulus),
% v (Poisson's ratio) and rho (mass density).
% 'Ke' is the sum of the flexural matrix and shear matrix.

% the order of the element's dof are like: [node_1_dof_1 ... node_1_dof_k ...... node_n_dof_1 ... node_n_dof_k],
% where n is the number of nodes and k is the number of dofs per node.

import FEM.*
syms x y;
E = material.E; v = material.v;
dx = element.dims.width;
dy = element.dims.height;
dz = element.dims.thickness;
N = getSF(element);
Cf = E/(1-v^2)*[1 v 0
                v 1 0
                0 0 (1-v)/2]*dz^3/12; % flexural constitutive matrix
Cs = 5/6*E/(2*(1+v))*[1 0
                      0 1]*dz; % shear constitutive matrix
switch element.type
    % Kirchhoff
    case {'ACM', 'BMF'}
        Bf = [-diff(diff(N, x), x)
              -diff(diff(N, y), y)
              -2*diff(diff(N, x), y)]; % strain-displacement matrix
        Ke = double(int(int(Bf'*Cf*Bf, y, -dy/2, dy/2), x, -dx/2, dx/2));
    % Mindlin
    case 'MB4'
        Bf = [diff(N(2,:), x)
              diff(N(3,:), y)
              diff(N(2,:), y) + diff(N(3,:), x)]; % flexural strain-displacement matrix
        Bs = [diff(N(1,:), x)
              diff(N(1,:), y)] + N(2:3,:); % shear strain-displacement matrix
        Ke = double(int(int(Bf'*Cf*Bf + Bs'*Cs*Bs, y, -dy/2, dy/2), x, -dx/2, dx/2));
end