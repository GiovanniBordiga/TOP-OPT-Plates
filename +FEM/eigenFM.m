function [eigenF, eigenM] = eigenFM(problem, element, x, PenK, PenM, optFindex, beta)
% Solve the eigenvalues problem K - omega^2*M = 0.
% 'problem' is a Problem object.
% 'nelx' and 'nely' are the number of elements along the two dimensions.
% 'element' is a FE object.
% 'x' is a nely-by-nelx matrix representing the density field on the plate.
% 'PenK' and 'PenM' are the penalization coefficients used in the SIMP model.
% 'nModes' is the number of eigenmodes (or eigenfrequencies) to return.
% 'eigenF' is the vector of the 'nModes' eigenvalues.
% 'eigenM' is the matrix of the 'nModes' eigenvectors (stored as columns).
% 'optFindex' is the eigenvalue's index to optimize.
% 'beta' is a stabilization factor for the mass matrix.

% The global numbering of the plate's dofs is ordered by columns.

nelx = problem.nelx;
nely = problem.nely;
n = element.nodes;     % nodes
ndof = element.ndof;   % dofs per node
Ke = element.K;
Me = element.M;
DOFindex = problem.DOFindex;

% assemble K and M assuming the same size for all the elements
x = reshape(x, 1, nelx*nely);
rowindex = kron(DOFindex, ones(1, ndof*n));
colindex = reshape(kron(reshape(DOFindex, ndof*n, nelx*nely), ones(1, ndof*n)), 1, nelx*nely*(ndof*n)^2);
K = sparse(rowindex, colindex, kron(x.^PenK, reshape(Ke, 1, (ndof*n)^2)));
M = sparse(rowindex, colindex, kron(x.^PenM, reshape(Me, 1, (ndof*n)^2)));
Mfake = sparse(rowindex, colindex, kron(x.^PenM, reshape(eye(ndof*n)/4, 1, (ndof*n)^2)));
M = M*(1-beta) + Mfake*beta;    % altering the mass matrix

% solve the eigenvalues problem
eigenM = zeros(ndof*(nely+1)*(nelx+1), optFindex);
freedof = problem.freedof;
[V, D] = eigs(K(freedof, freedof), M(freedof, freedof), optFindex, 'sm'); % returned the 'nModes' smallest eigenvalues and relative eigenvectors (normalized to unit modal masses)
eigenF = diag(D);
eigenM(freedof,:) = V;
% return sorted eigenvectors and eigenvalues
for i = 1:length(eigenF)
    for j = i:length(eigenF)
        if eigenF(j) < eigenF(i)
            temp = eigenF(i);
            eigenF(i) = eigenF(j);
            eigenF(j) = temp;
            
            temp = eigenM(:, i);
            eigenM(:, i) = eigenM(:, j);
            eigenM(:, j) = temp;
        end
    end
end
end

