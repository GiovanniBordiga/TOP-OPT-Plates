function [Kf, Ks] = getK(element, material)
% Return the stiffness matrices of the given element.
% 'element' is a FE object.
% 'material' is a struct containing two attributes: E (Young's modulus) and
% v (Poisson's ratio).
% 'Kf' is the matrix of the flexural elastic energy.
% 'Ks' is the matrix of the shear elastic energy (zero for Kirchhoff elements).

% the order of the element's dof are like: [node_1_dof_1 ... node_1_dof_k ...... node_n_dof_1 ... node_n_dof_k],
% where n is the number of nodes and k is the number of dofs per node.

import FEM.*
syms x y z;
E = material.E; v = material.v;
dims = element.getDims();
N = getSF(element);
Cf = E/(1-v^2)*[1 v 0
                v 1 0
                0 0 (1-v)/2]*dims.thickness^3/12; % flexular constitutive matrix
Cs = 5/6*E/(2*(1+v))*[1 0
                      0 1]*dims.thickness; % shear constitutive matrix
switch element.getType()
    % Kirchhoff
    case {'ACM', 'BMF'}
        Bf = [-diff(diff(N,x),x)
             -diff(diff(N,y),y)
             -2*diff(diff(N,x),y)]; % strain-displacement matrix
        Kf = double(int(int(Bf'*Cf*Bf,y,-dims.height/2,dims.height/2),x,-dims.width/2,dims.width/2));
        Ks = 0;
    % Mindlin
    case 'MB4'
        Bf = [diff(N(2,:),x)
              diff(N(3,:),x)
              diff(N(2,:),y) + diff(N(3,:),x)]; % flexural strain-displacement matrix
        Bs = [diff(N(1,:),x)
              diff(N(1,:),y)] + N(2:3,:); % shear strain-displacement matrix
        Kf = double(int(int(Bf'*Cf*Bf,y,-dims.height/2,dims.height/2),x,-dims.width/2,dims.width/2));
        Ks = double(int(int(Bs'*Cs*Bs,y,-dims.height/2,dims.height/2),x,-dims.width/2,dims.width/2));
end