function U = FEM(problem, nelx, nely, element, material, x, CoPen)
% Solve the system KU=F.
% 'nelx' and 'nely' are the number of elements along the two dimensions.
% 'element' is a FE object.
% 'material' is a struct containing two attributes: E (Young's modulus) and
% v (Poisson's ratio).
% 'x' is a nely-by-nelx matrix representing the density field on the plate.
% 'U' is the global dofs vector.

% the order of the element's dof are like: [node_1_dof_1 ... node_1_dof_k ...... node_n_dof_1 ... node_n_dof_k],
% where n is the number of nodes and k is the number of dofs per node.

import FEM.*
n = element.getNodes();     % nodes
ndof = element.getNDof();   % dofs per node
F = sparse(ndof*(nely+1)*(nelx+1), 1);  
U = zeros(ndof*(nely+1)*(nelx+1), 1);

[Kf, Ks] = getK(element, material);
% assemble K assuming the same size for all the elements
% probably not so efficient (TEST!)
% for elx = 1:nelx
%     for ely = 1:nely
%         nodesnum = [(elx-1)*(nely+1) + ely + 1
%                     (elx)*(nely+1) + ely + 1
%                     (elx)*(nely+1) + ely
%                     (elx-1)*(nely+1) + ely];
%         % create global dof index for the element
%         dofindex = [];
%         for i = 1:n
%             dofindex = union(dofindex, nodesnum(i):nodesnum(i)+ndof-1);
%         end
%         K(dofindex, dofindex)= K(dofindex, dofindex) + x(ely, elx)*(Kf + Ks); 
%     end
% end

% create global dof index
DOFindex = [];
for elx = 1:nelx
    for ely = 1:nely
        nodesnum = [(elx-1)*(nely+1) + ely + 1
                    (elx)*(nely+1) + ely + 1
                    (elx)*(nely+1) + ely
                    (elx-1)*(nely+1) + ely];    % global nodes numbers of the current element
        for i = 1:n
            DOFindex = cat(2, DOFindex, (nodesnum(i)-1)*ndof+1:nodesnum(i)*ndof);
        end
    end
end
% assemble K assuming the same size for all the elements
x = reshape(x, 1, nelx*nely);
rowindex = kron(DOFindex, ones(1, ndof*n));
colindex = reshape(kron(reshape(DOFindex, ndof*n, nelx*nely), ones(1, ndof*n)), 1, nelx*nely*(ndof*n)^2);
K = sparse(rowindex, colindex, kron(x.^CoPen, reshape(Kf+Ks, 1, (ndof*n)^2)));

% loads and constraints
switch problem
    case 'test'
        F(1:ndof:ndof*(nely+1)*(nelx+1)) = 1;   % distributed load
        alldof = 1:ndof*(nelx+1)*(nely+1);
        fixeddof = 1:ndof*(nely+1);             % one edge clamped
        freedof = setdiff(alldof, fixeddof);
        U(fixeddof) = 0;
    case 'a'
        % <https://goo.gl/FQdEKB link to image of the problem>
        F(ndof*(nely+1)*nelx/2 + 1) = 1;        % concentrated load
        alldof = 1:ndof*(nelx+1)*(nely+1);
        left = 1:ndof:ndof*nely+1;              % supported left edge
        bottom = ndof*nely+1:ndof*(nely+1):ndof*((nely+1)*(nelx+1)-1)+1; % supported bottom edge
        right = ndof*(nely+1)*nelx+1:ndof:ndof*((nely+1)*(nelx+1)-1)+1; % supported right edge
        fixeddof = union(union(left, bottom), right);
        freedof = setdiff(alldof, fixeddof);    % free top edge
        U(fixeddof) = 0;
end

% solve the system
U(freedof) = K(freedof, freedof) \ F(freedof);
end