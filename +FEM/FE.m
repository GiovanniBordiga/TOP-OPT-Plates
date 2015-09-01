classdef FE
    %FE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
        dims
    end
    
    methods
        function obj = FE(type, dims)
            obj.type = type;
            obj.dims = dims;
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
    end
    
end

