classdef Problem
    % Class representing a FEM problem.
    % Given the plate's dimensions, a finite element and a 'problemId', it
    % defines the boundary conditions in its properties: 'fixeddof',
    % 'freedof' and 'F'.
    
    properties
        nelx            % number of elements along the x axis
        nely            % number of elements along the y axis
        problemId       % string identifing certain boundary conditions
        fixeddof        % idexes of the constrained dofs
        freedof         % idexes of the free dofs
        F               % global vector of loads
    end
    
    methods
        function obj = Problem(nelx, nely, element, problemId)
            obj.nelx = nelx;
            obj.nely = nely;
            obj.problemId = problemId;
            ndof = element.ndof;
            % loads and constraints
            obj.F = sparse(ndof*(nely+1)*(nelx+1), 1);
            switch obj.problemId
                case 'test1'
                    obj.F(1:ndof:ndof*(nely+1)*(nelx+1)) = 0.001;   % ~distributed load
                    alldof = 1:ndof*(nelx+1)*(nely+1);
                    obj.fixeddof = 1:ndof*(nely+1);             % one edge clamped
                case 'test2'
                    obj.F(ndof*((nely+1)*nelx/2 + nely/2) + 1) = 1; % centered concentrated load
                    alldof = 1:ndof*(nelx+1)*(nely+1);
                    left = 1:ndof:ndof*nely+1;              % supported left edge
                    right = ndof*(nely+1)*nelx+1:ndof:ndof*((nely+1)*(nelx+1)-1)+1; % supported right edge
                    obj.fixeddof = union(left, right);
                case 'test3'
                    obj.F(ndof*((nely+1)*nelx/2 + nely/2) + 1) = 1; % centered concentrated load
                    alldof = 1:ndof*(nelx+1)*(nely+1);
                    left = 1:ndof:ndof*nely+1;              % supported left edge
                    bottom = ndof*nely+1:ndof*(nely+1):ndof*((nely+1)*(nelx+1)-1)+1; % supported bottom edge
                    right = ndof*(nely+1)*nelx+1:ndof:ndof*((nely+1)*(nelx+1)-1)+1; % supported right edge
                    top = 1:ndof*(nely+1):ndof*(nely+1)*nelx+1; % supported top edge
                    obj.fixeddof = union(union(union(left, right), bottom), top);
                case 'a'
                    % <https://goo.gl/FQdEKB link to image of the problem>
                    obj.F(ndof*(nely+1)*nelx/2 + 1) = 1;        % concentrated load
                    alldof = 1:ndof*(nelx+1)*(nely+1);
                    left = 1:ndof:ndof*nely+1;              % supported left edge
                    bottom = ndof*nely+1:ndof*(nely+1):ndof*((nely+1)*(nelx+1)-1)+1; % supported bottom edge
                    right = ndof*(nely+1)*nelx+1:ndof:ndof*((nely+1)*(nelx+1)-1)+1; % supported right edge
                    obj.fixeddof = union(union(left, bottom), right);
            end
            obj.freedof = setdiff(alldof, obj.fixeddof);
        end
    end
    
end

