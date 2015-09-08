function plotDeformed(nelx, nely, element, xpen, Ke, U)
% Plot the deformed configuration of the plate colored with the
% distribution of the elastic energy.
% 'nelx' and 'nely' are the number of elements along the two dimensions.
% 'element' is a FE object.
% 'U' is the global dof vector returned by FEM function.

faces = [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8]; % element's faces for patch
n = element.getNodes();     % nodes
ndof = element.getNDof();   % dofs per node
dx = element.getWidth();
dy = element.getHeight();
dz = element.getThickness();
figure(4);
colorbar;
axis equal; axis tight;         % set axis to match plate dimensions
title('Deformed configuration and elastic energy distribution'); % set title
xlabel('x'); ylabel('y');       % set axis labels
view(3);                        % set 3D view
hold on;
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
        xe = [(elx-1)*dx elx*dx elx*dx (elx-1)*dx (elx-1)*dx elx*dx elx*dx (elx-1)*dx]';
        ye = [ely*dy ely*dy (ely-1)*dy (ely-1)*dy ely*dy ely*dy (ely-1)*dy (ely-1)*dy]';
        ze = [ones(1,4)*(-dz/2) ones(1,4)*dz/2]';
        undef = [xe ye ze]; % undeformed vertices coordinates
        disp = [[Ue(2:ndof:end); Ue(2:ndof:end)].*ze [Ue(3:ndof:end); Ue(3:ndof:end)].*ze [Ue(1:ndof:end); Ue(1:ndof:end)]]; % vertices displacements
        def = undef + disp; % deformed vertices coordinates
        patch('Faces', faces, 'Vertices', def, 'FaceColor', 'flat',...
            'CDataMapping', 'scaled', 'EdgeColor', [0.1 0.1 0.1],...
            'FaceVertexCData', 0.5*xpen(ely, elx)*Ue'*Ke*Ue);
    end
end
hold off;
end