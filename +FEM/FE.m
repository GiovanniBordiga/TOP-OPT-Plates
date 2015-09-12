classdef FE
    % Class representing a finite element.
    % This class supports three types of finite element, two Kirchhoff
    % elements (ACM and BMF) and one Mindlin element (MB4).
    % The order of the element's dof are like: [node_1_dof_1 ... node_1_dof_k ...... node_n_dof_1 ... node_n_dof_k],
    % where n is the number of nodes and k is the number of dofs per node.
    % The numeration of the element's nodes is anticlockwise starting from
    % the bottom left corner.
    
    properties
        type        % string representing the element type.
        dims        % struct which defines the element's dimensions: width, height, thickness.
        nodes       % number of nodes.
        ndof        % number of dofs per node.
        material    % struct containing three attributes: E (Young's modulus),
                    % v (Poisson's ratio) and rho (mass density).
        N           % shape functions m-by-n matrix, where n is the number
                    % of fields modeled and n is the number of dofs per element.
        K           % stiffness matrix.
        M           % mass matrix.
    end
    
    methods
        % Constructor: call this method to create a finite element object
        function element = FE(type, dims, material)
            element.type = type;
            element.dims = dims;
            element.material = material;
            switch type
                % Kirchhoff
                case 'ACM'
                    element.nodes = 4;
                    element.ndof = 3;
                case 'BMF'
                    element.nodes = 4;
                    element.ndof = 4;
                % Mindlin
                case 'MB4'      % Mindlin bilinear 4 nodes
                    element.nodes = 4;
                    element.ndof = 3;
            end
            element.N = element.buildSF();
            element.K = element.buildK();
            element.M = element.buildM();
        end
        
        % Shape functions builder
        function N = buildSF(element)
            syms x y;
            dims = element.dims;   % element's dimensions
            n = element.nodes;     % nodes
            ndof = element.ndof;   % dofs per node
            coord = [-dims.width/2 -dims.height/2
                      dims.width/2 -dims.height/2
                      dims.width/2  dims.height/2
                     -dims.width/2  dims.height/2]; % nodes coordinates
            switch element.type
                % Kirchhoff
                case 'ACM'
                    P = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3]'; % w=P'a
                    Px = diff(P,x); Py = diff(P,y);
                    A = zeros(n*ndof); % Aa=w
                    for i = 1:n
                        A(ndof*(i-1)+1:ndof*i,:) = subs([P'; Px'; Py'], {x,y}, {coord(i,1), coord(i,2)});
                    end
                    N = (P'*A^-1); % Nw=P'(A^-1)w
                case 'BMF'
                    P = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x^2*y^2 x*y^3 x^3*y^2 x^2*y^3 x^3*y^3]';
                    Px = diff(P,x); Py = diff(P,y); Pxy = diff(Px,y);
                    A = zeros(n*ndof); % Aa=w
                    for i = 1:n
                        A(ndof*(i-1)+1:ndof*i,:) = subs([P'; Px'; Py'; Pxy'], {x,y}, {coord(i,1), coord(i,2)});
                    end
                    N = (P'*A^-1); % Nw=P'(A^-1)w
                % Mindlin
                case 'MB4'          % Mindlin bilinear 4 nodes
                    P = [1 x y x*y]';
                    A = zeros(n); % Aa=w
                    for i = 1:n
                        A(i,:) = subs(P, {x,y}, {coord(i,1), coord(i,2)})';
                    end
                    psi = (P'*A^-1); % Nw=P'(A^-1)w
                    for i = 1:n
                        N(:, ndof*(i-1)+1:ndof*i) = eye(ndof)*psi(i);
                    end
            end
        end
        
        % Stiffness matrix builder
        function K = buildK(element)
            syms x y;
            E = element.material.E; v = element.material.v;
            dx = element.dims.width;
            dy = element.dims.height;
            dz = element.dims.thickness;
            N = element.N;
            Cf = E/(1-v^2)*[1 v 0
                            v 1 0
                            0 0 (1-v)/2]*dz^3/12; % flexural constitutive matrix
            Cs = 5/6*E/(2*(1+v))*[1 0
                                  0 1]*dz; % shear constitutive matrix
            switch element.type
                % Kirchhoff
                case {'ACM', 'BMF'}
                    Bf = [-diff(diff(N, x), x)
                          -diff(diff(N, y), y)
                          -2*diff(diff(N, x), y)]; % strain-displacement matrix
                    K = double(int(int(Bf'*Cf*Bf, y, -dy/2, dy/2), x, -dx/2, dx/2));
                % Mindlin
                case 'MB4'
                    Bf = [diff(N(2,:), x)
                          diff(N(3,:), y)
                          diff(N(2,:), y) + diff(N(3,:), x)]; % flexural strain-displacement matrix
                    Bs = [diff(N(1,:), x)
                          diff(N(1,:), y)] + N(2:3,:); % shear strain-displacement matrix
                    K = double(int(int(Bf'*Cf*Bf + Bs'*Cs*Bs, y, -dy/2, dy/2), x, -dx/2, dx/2));
            end
        end
        
        % Mass matrix builder
        function M = buildM(element)
            syms x y;
            dx = element.dims.width;
            dy = element.dims.height;
            dz = element.dims.thickness;
            rho = element.material.rho;         % mass density
            N = element.N;
            C = [dz 0 0
                0 dz^3/12 0
                0 0 dz^3/12]*rho;
            switch element.type
                % Kirchhoff
                case {'ACM', 'BMF'}
                    B = [N; -diff(N, x); -diff(N, y)]; %#ok<*PROP>
                    M = double(int(int(B'*C*B, y, -dy/2, dy/2), x, -dx/2, dx/2));
                % Mindlin
                case 'MB4'
                    M = double(int(int(N'*C*N, y, -dy/2, dy/2), x, -dx/2, dx/2));
            end
        end
    end
    
end