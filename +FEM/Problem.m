classdef Problem
    % Class representing a FEM problem.
    % Given the plate's dimensions, a finite element and a 'problemId', it
    % defines the boundary conditions in its properties: 'fixeddof',
    % 'freedof' and 'F'.
    
    properties
        nelx            % number of elements along the x axis
        nely            % number of elements along the y axis
        problemId       % string identifing certain boundary conditions
        DOFindex        % global dofs' index
        fixeddof        % idexes of the constrained dofs
        freedof         % idexes of the free dofs
        F               % global vector of loads
    end
    
    methods
        function obj = Problem(nelx, nely, element, problemId)
            obj.nelx = nelx;
            obj.nely = nely;
            obj.problemId = problemId;
            n = element.nodes;
            ndof = element.ndof;
            
            % loads and constraints
            obj.F = sparse(ndof*(nely+1)*(nelx+1), 1);
            alldof = 1:ndof*(nelx+1)*(nely+1);
            switch obj.problemId
            % <https://goo.gl/v4SDZr link to image of the problems>
                case 'test1'
                    obj.F(1:ndof:ndof*(nely+1)*(nelx+1)) = 0.001;   % ~distributed load
                    obj.fixeddof = 1:ndof*(nely+1);                 % left edge clamped
                case 'test2'
                    obj.F(ndof*((nely+1)*nelx/2 + nely/2) + 1) = 1; % centered concentrated load
                    left = 1:ndof:ndof*nely+1;                      % supported left edge
                    right = ndof*(nely+1)*nelx+1:ndof:ndof*((nely+1)*(nelx+1)-1)+1; % supported right edge
                    obj.fixeddof = union(left, right);
                case 'test3'
                    obj.F(ndof*((nely+1)*nelx/2 + nely/2) + 1) = 1; % centered concentrated load
                    left = 1:ndof:ndof*nely+1;                      % supported left edge
                    bottom = ndof*nely+1:ndof*(nely+1):ndof*((nely+1)*(nelx+1)-1)+1; % supported bottom edge
                    right = ndof*(nely+1)*nelx+1:ndof:ndof*((nely+1)*(nelx+1)-1)+1;  % supported right edge
                    top = 1:ndof*(nely+1):ndof*(nely+1)*nelx+1;     % supported top edge
                    obj.fixeddof = union(union(union(left, right), bottom), top);
                case 'a'
                    obj.F(ndof*(nely+1)*nelx/2 + 1) = 1;            % concentrated load
                    left = 1:ndof:ndof*nely+1;                      % supported left edge
                    bottom = ndof*nely+1:ndof*(nely+1):ndof*((nely+1)*(nelx+1)-1)+1; % supported bottom edge
                    right = ndof*(nely+1)*nelx+1:ndof:ndof*((nely+1)*(nelx+1)-1)+1;  % supported right edge
                    obj.fixeddof = union(union(left, bottom), right);
                case 'b'
                    % same as case 'a' with clamped edge
                    obj.F(ndof*(nely+1)*nelx/2 + 1) = 1;        % concentrated load
                    left = 1:ndof*(nely+1);                       % clamped left edge
                    bottom = [];
                    for i = 1:ndof
                        bottom_new = ndof*nely+i:ndof*(nely+1):ndof*((nely+1)*(nelx+1)-1)+i;     % clamped bottom edge
                        bottom = union(bottom_new, bottom);
                    end
                    right = ndof*(nely+1)*nelx+1:ndof*(nely+1)*(nelx+1);                   % supported right edge
                    obj.fixeddof = union(union(left, bottom), right);
                case 'c'
                    obj.F(ndof*(nely+1)*nelx/2 + 1) = 1;            % concentrated load in the middle
                    left = 1:ndof:ndof*nely+1;                      % supported left edge
                    right = ndof*(nely+1)*nelx+1:ndof:ndof*((nely+1)*(nelx+1)-1)+1; % supported right edge
                    obj.fixeddof = union(left, right);
                case 'd'
                    obj.F(1) = 1;                                   % concentrated load on the left corner
                    obj.F((ndof*(nely+1)*nelx/2) + 1) = 1;          % concentrated load in the middle
                    obj.F((ndof*(nely+1)*nelx) + 1) = 1;            % concentrated load on right corner
                    bottom = [];
                    for i = 1:ndof
                        bottom_new = ndof*nely+i:ndof*(nely+1):ndof*((nely+1)*(nelx+1)-1)+i;     % clamped bottom edge
                        bottom = union(bottom_new, bottom);
                    end
                    obj.fixeddof = bottom;
                case 'e'
                    obj.F(ndof*(nely+1)*nelx + 1) = 1;              % concentrated load right corner
                    obj.fixeddof = 1:ndof*(nely+1);                 % left edge clamped
            end
            obj.freedof = setdiff(alldof, obj.fixeddof);
            
            % create global dofs' index
            obj.DOFindex = [];
            for elx = 1:nelx
                for ely = 1:nely
                    nodesnum = [(elx-1)*(nely+1) + ely + 1
                                (elx)*(nely+1) + ely + 1
                                (elx)*(nely+1) + ely
                                (elx-1)*(nely+1) + ely];    % global nodes numbers of the current element
                    for i = 1:n
                        obj.DOFindex = cat(2, obj.DOFindex, (nodesnum(i)-1)*ndof+1:nodesnum(i)*ndof);
                    end
                end
            end
        end
    end
    
end

