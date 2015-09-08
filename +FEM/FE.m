classdef FE
    % Class representing a finite element.
    
    properties
        type        % string representing the element type
        dims        % struct which defines the element's dimensions: width, height, thickness
        nodes       % number of nodes
        ndof        % number of dofs per node
    end
    
    methods
        function obj = FE(type, dims)
            obj.type = type;
            obj.dims = dims;
            switch type
                % Kirchhoff
                case 'ACM'
                    obj.nodes = 4;
                    obj.ndof = 3;
                case 'BMF'
                    obj.nodes = 4;
                    obj.ndof = 4;
                % Mindlin
                case 'MB4'      % Mindlin bilinear 4 nodes
                    obj.nodes = 4;
                    obj.ndof = 3;
            end
        end
    end
    
end

