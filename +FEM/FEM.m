function U = FEM(problem, nelx, nely, element, x, CoPen)
% Solve the system KU=F.
% 'problem' is a Problem object.
% 'nelx' and 'nely' are the number of elements along the two dimensions.
% 'element' is a FE object.
% 'x' is a nely-by-nelx matrix representing the density field on the plate.
% 'CoPen' is the penalization coefficient used in the SIMP model.
% 'U' is the global dofs vector.

% The global numbering of the plate's dofs is ordered by columns.

n = element.nodes;     % nodes
ndof = element.ndof;   % dofs per node
Ke = element.K;

U = zeros(ndof*(nely+1)*(nelx+1), 1);

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
K = sparse(rowindex, colindex, kron(x.^CoPen, reshape(Ke, 1, (ndof*n)^2)));

% solve the system
freedof = problem.freedof;
F = problem.F;
U(freedof) = K(freedof, freedof) \ F(freedof);
end