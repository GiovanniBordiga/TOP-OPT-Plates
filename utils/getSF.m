function N = getSF(element, dims)
syms x y;
switch element
    % Kirchhoff
    case 'ACM'
        n = 4;              % nodes
        ndof = 3;           % dof per node
        P = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3]'; % w=P'a
        Px = diff(P,x); Py = diff(P,y);
        coord = [-dims.width/2 -dims.height/2
                  dims.width/2 -dims.height/2
                  dims.width/2  dims.height/2
                 -dims.width/2  dims.height/2]; % nodes coordinates
        A = zeros(n*ndof); % Aa=w
        for i = 1:n
            A(ndof*(i-1)+1:ndof*i,:) = subs([P'; Px'; Py'], {x,y}, {coord(i,1), coord(i,2)});
        end
        N=(P'*A^-1); % Nw=P'(A^-1)w
    case 'BMF'
        n = 4;
        ndof = 4;
        P = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x^2*y^2 x*y^3 x^3*y^2 x^2*y^3 x^3*y^3]';
        Px = diff(P,x); Py = diff(P,y); Pxy = diff(Px,y);
        coord = [-dims.width/2 -dims.height/2
                  dims.width/2 -dims.height/2
                  dims.width/2  dims.height/2
                 -dims.width/2  dims.height/2]; % nodes coordinates
        A = zeros(n*ndof); % Aa=w
        for i = 1:n
            A(ndof*(i-1)+1:ndof*i,:) = subs([P'; Px'; Py'; Pxy'], {x,y}, {coord(i,1), coord(i,2)});
        end
        N=(P'*A^-1); % Nw=P'(A^-1)w
    % Mindlin
    case 'MB4'          % Mindlin bilinear 4 nodes
        n = 4;
        ndof = 3;
        P = [1 x y x*y]';
        coord = [-dims.width/2 -dims.height/2
                  dims.width/2 -dims.height/2
                  dims.width/2  dims.height/2
                 -dims.width/2  dims.height/2]; % nodes coordinates
        A = zeros(n); % Aa=w
        for i = 1:n
            A(i,:) = subs(P, {x,y}, {coord(i,1), coord(i,2)})';
        end
        psi=(P'*A^-1); % Nw=P'(A^-1)w
        for i = 1:n
            N(:, ndof*(i-1)+1:ndof*i) = eye(ndof)*psi(i);
        end
end
end