function U = FEM(nelx, nely, element, dims, material, x)
ndof = 3; % TODO to get from element?
n = 4;    % TODO to get from element?
F = sparse(ndof*(nely+1)*(nelx+1), 1);  
U = zeros(ndof*(nely+1)*(nelx+1), 1);

[Kf, Ks] = getK(element, dims, material);
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
                    (elx-1)*(nely+1) + ely];
        for i = 1:n
            DOFindex = union(DOFindex, nodesnum(i):nodesnum(i)+ndof-1);
        end
    end
end
% assemble K assuming the same size for all the elements
x = reshape(x', 1, nelx*nely);
rowindex = kron(DOFindex, ones(1, ndof*n));
colindex = reshape(kron(reshape(DOFindex, ndof*n, nelx*nely), ones(1, ndof*n)), 1, nelx*nely*(ndof*n)^2);
K = sparse(rowindex, colindex, kron(x, reshape(Kf+Ks, 1, (ndof*n)^2)));

% define loads and constraints - MODIFY AS YOU LIKE
F(1:3:3*(nely+1)*(nelx+1),1) = 1;
alldof = 1:ndof*(nelx+1)*(nely+1);
fixeddof = 1:ndof*(nely+1); % one edge clamped
freedof = setdiff(alldof, fixeddof);
U(fixeddof) = 0;

% solve the system
U(freedof) = K(freedof, freedof) \ F(freedof);
end