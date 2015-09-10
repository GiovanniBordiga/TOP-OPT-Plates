function [eigenF, eigenM] = eigenFM(problem, nelx, nely, element, material, x, PenK, PenM, nModes)
% Solve the eigenvalues problem K - omega^2*M = 0.
% 'problem' is a Problem object.
% 'nelx' and 'nely' are the number of elements along the two dimensions.
% 'element' is a FE object.
% 'material' is a struct containing three attributes: E (Young's modulus),
% v (Poisson's ratio) and rho (mass density).
% 'x' is a nely-by-nelx matrix representing the density field on the plate.
% 'PenK' and 'PenM' are the penalization coefficients used in the SIMP model.
% 'nModes' is the number of eigenmodes (or eigenfrequencies) to return.
% 'eigenF' is the vector of the 'nModes' eigenvalues.
% 'eigenM' is the matrix of the 'nModes' eigenvectors (stored as columns).

import FEM.*

n = element.nodes;     % nodes
ndof = element.ndof;   % dofs per node
Ke = getK(element, material);
Me = getM(element, material);

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
% assemble K and M assuming the same size for all the elements
x = reshape(x, 1, nelx*nely);
% PenM = x; PenM(PenM<0.1) = 6; PenM(PenM~=6) = 3;
rowindex = kron(DOFindex, ones(1, ndof*n));
colindex = reshape(kron(reshape(DOFindex, ndof*n, nelx*nely), ones(1, ndof*n)), 1, nelx*nely*(ndof*n)^2);
K = sparse(rowindex, colindex, kron(x.^PenK, reshape(Ke, 1, (ndof*n)^2)));
M = sparse(rowindex, colindex, kron(x.^PenM, reshape(Me, 1, (ndof*n)^2)));

% solve the eigenvalues problem
eigenM = zeros(ndof*(nely+1)*(nelx+1), nModes);
freedof = problem.freedof;
[V, D] = eigs(K(freedof, freedof), M(freedof, freedof), nModes, 'sm'); % returned the 'nModes' smallest eigenvalues and relative eigenvectors (normalized to unit modulus)
eigenF = diag(D);
eigenM(freedof,:) = V;
end

