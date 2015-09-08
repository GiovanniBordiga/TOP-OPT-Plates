function M = getM(element, material)
% Return the mass matrix of the given element (considering both the
% traslational and rotational inertia).
% 'element' is a FE object.
% 'material' is a struct containing three attributes: E (Young's modulus),
% v (Poisson's ratio) and rho (mass density).

import FEM.*
syms x y;
dx = element.dims.width;
dy = element.dims.height;
dz = element.dims.thickness;
rho = material.rho;         % mass density
N = getSF(element);
C = [dz 0 0
     0 dz^3/12 0
     0 0 dz^3/12]*rho;
switch element.getType()
    % Kirchhoff
    case {'ACM', 'BMF'}
        B = [N; -diff(N, x); -diff(N, y)];
        M = double(int(int(B'*C*B, y, -dy/2, dy/2), x, -dx/2, dx/2));
    % Mindlin
    case 'MB4'
        M = double(int(int(N'*C*N, y, -dy/2, dy/2), x, -dx/2, dx/2));
end
end

