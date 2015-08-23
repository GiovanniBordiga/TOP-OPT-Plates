function U = FEM(nelx, nely, element, dims, material, x)
ngdl = 3; % TODO to get from element?
n = 4;    % TODO to get from element?
K = sparse(ngdl*(nelx+1)*(nely+1), ngdl*(nelx+1)*(nely+1));
F = sparse(ngdl*(nely+1)*(nelx+1), 1);  
U = zeros(ngdl*(nely+1)*(nelx+1), 1);

[Kf, Ks] = getK(element, dims, material);
% assemble K assuming the same size for all the elements
% TODO better to optimize using sparse!
for elx = 1:nelx
    for ely = 1:nely
        nodesnum = [(elx-1)*(nely+1) + ely + 1
                    (elx)*(nely+1) + ely + 1
                    (elx)*(nely+1) + ely
                    (elx-1)*(nely+1) + ely];
        % create global dof index
        for i = 1:n
            index = union(index, nodesnum(i):nodesnum(i)+ngdl);
        end
        K(index, index)= K(index, index) + x(ely, elx)*(Kf + Ks); 
    end
end

% define loads and constraints - MODIFY AS YOU LIKE
F(1:3:3*(nely+1)*(nelx+1),1) = 1;
alldof = 1:ngdl*(nelx+1)*(nely+1);
fixeddof = 1:ngdl*(nely+1); % one edge clamped
freedof = setdiff(alldof, fixeddof);
U(fixeddof) = 0;

% solve the system
U(freedof) = K(freedof, freedof) \ F(freedof);
end