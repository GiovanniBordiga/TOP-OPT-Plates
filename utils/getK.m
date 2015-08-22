function [Kf, Ks] = getK(element, dims, E, v)
syms x y z;
N = getSF(element,dims);
switch element
    % Kirchhoff
    case {'ACM', 'BMF'}
        C = E/(1-v^2)*[1 v 0
            v 1 0
            0 0 (1-v)/2]*dims.thickness^3/12;
        B = [-diff(diff(N,x),x)'
             -diff(diff(N,y),y)'
             -2*diff(diff(N,x),y)']; % strain-displacement matrix
        Kf = double(int(int(B'*C*B,y,-dims.height/2,dims.height/2),x,-dims.width/2,dims.width/2));
        Ks = 0;
    % Mindlin
    case 'MB4'
        % TODO
end