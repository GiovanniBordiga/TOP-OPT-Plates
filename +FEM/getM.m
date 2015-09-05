function M = getM(element, material)
% Return the mass matrix of the given element (considering both the
% traslational and rotational inertia).
% 'element' is a FE object.
% 'material' is a struct containing three attributes: E (Young's modulus),
% v (Poisson's ratio) and rho (mass density).

import FEM.*
syms x y;
dx = element.getWidth();
dy = element.getHeight();
dz = element.getThickness();
rho = material.rho;         % mass density
N = getSF(element);
switch element.getType()
    case {'ACM', 'BMF'}
        % TODO thinking...
    case 'MB4'
        C = [dz 0 0
             0 dz^3/12 0
             0 0 dz^3/12]*rho;
        M = double(int(int(N'*C*N, y, -dy/2, dy/2), x, -dx/2, dx/2));
end
end

