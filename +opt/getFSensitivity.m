function dF = getFSensitivity(nelx, nely, element, x, PenK, PenM, Ke, Me, eigenF, eigenM)
% Return the sensitivity of the lowest eigenfrequency with respect to densities.
% 'nelx' and 'nely' are the number of elements along the two dimensions.
% 'element' is a FE object.
% 'x' is a nely-by-nelx matrix representing the density field on the plate.
% 'PenK' and 'PenM' are the penalization coefficients used in the SIMP model.
% 'Ke' and 'Me' are, respectively, the stiffness matrix and the mass matrix
% of an element (assuming that is the same for every element).
% 'eigenF' is the vector of the eigenvalues.
% 'eigenM' is the matrix of the eigenvectors (stored as columns).

n = element.nodes;     % nodes
ndof = element.ndof;   % dofs per node
dF = zeros(nely, nelx);
lowF = min(eigenF);         % get the lowest eigenfrequency
lowM = eigenM(:, eigenF==lowF);      % get the corresponding eigenmode
for elx = 1:nelx
    for ely = 1:nely
        eDOFindex = [];     % global dofs index of the current element, needed to access the element's dofs from 'lowM'
        nodesnum = [(elx-1)*(nely+1) + ely + 1
                    (elx)*(nely+1) + ely + 1
                    (elx)*(nely+1) + ely
                    (elx-1)*(nely+1) + ely];    % global nodes numbers of the current element
        for i = 1:n
            eDOFindex = cat(2, eDOFindex, (nodesnum(i)-1)*ndof+1:nodesnum(i)*ndof);
        end
        lowMe = lowM(eDOFindex);      % dofs of the current element
        dF(ely, elx) = lowMe'*(PenK*x(ely, elx)^(PenK-1)*Ke - lowF * PenM*x(ely, elx)^(PenM-1)*Me)*lowMe;
    end
end
end