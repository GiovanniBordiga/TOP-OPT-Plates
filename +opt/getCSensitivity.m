function [dC, C] = getCSensitivity(nelx, nely, element, x, CoPen, U)
% Return the sensitivity of the compliance with respect to densities.
% 'nelx' and 'nely' are the number of elements along the two dimensions.
% 'element' is a FE object.
% 'x' is a nely-by-nelx matrix representing the density field on the plate.
% 'CoPen' is the penalization coefficient used in the SIMP model.
% 'Ke' is the stiffness matrix of an element (assuming that is the same for
% every element).
% 'U' is the global dofs vector returned by FEM function.
% 'dC' and 'C' are, respectively, the sensitivity (stored as a matrix) and
% the total compliance.

n = element.nodes;     % nodes
ndof = element.ndof;   % dofs per node
Ke = element.K;
dC = zeros(nely, nelx);
C = 0;
for elx = 1:nelx
    for ely = 1:nely
        eDOFindex = [];     % global dofs index of the current element, needed to access the element's dofs from 'U'
        nodesnum = [(elx-1)*(nely+1) + ely + 1
                    (elx)*(nely+1) + ely + 1
                    (elx)*(nely+1) + ely
                    (elx-1)*(nely+1) + ely];    % global nodes numbers of the current element
        for i = 1:n
            eDOFindex = cat(2, eDOFindex, (nodesnum(i)-1)*ndof+1:nodesnum(i)*ndof);
        end
        Ue = U(eDOFindex);      % dofs of the current element
        dC(ely, elx) = - CoPen * x(ely, elx)^(CoPen-1) * Ue'*Ke*Ue;
        C = C + x(ely, elx)^CoPen * Ue'*Ke*Ue;
    end
end
end

