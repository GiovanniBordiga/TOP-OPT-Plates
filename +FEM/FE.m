classdef FE
    % Class representing a finite element.
    
    properties
        type            % string representing the element type
        dims            % struct which defines the element's dimensions: width, height, thickness
        nodes           % number of nodes
        ndof            % number of dof per node
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
                case 'MB4'
                    obj.nodes = 4;
                    obj.ndof = 3;
            end
        end
        function type = getType(obj)
            type = obj.type;
        end
        function dims = getDims(obj)
            dims = obj.dims;
        end
        function width = getWidth(obj)
            width = obj.dims.width;
        end
        function height = getHeight(obj)
            height = obj.dims.height;
        end
        function thickness = getThickness(obj)
            thickness = obj.dims.thickness;
        end
        function nodes = getNodes(obj)
            nodes = obj.nodes;
        end
        function ndof = getNDof(obj)
            ndof = obj.ndof;
        end
    end
    
end

