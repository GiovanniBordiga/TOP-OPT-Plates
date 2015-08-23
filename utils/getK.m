function [Kf, Ks] = getK(element, dims, E, v)
syms x y z;
N = getSF(element,dims);
Cf = E/(1-v^2)*[1 v 0
                v 1 0
                0 0 (1-v)/2]*dims.thickness^3/12; % flexular constitutive matrix
Cs = 5/6*E/(2*(1+v))*[1 0
                      0 1]*dims.thickness; % shear constitutive matrix
switch element
    % Kirchhoff
    case {'ACM', 'BMF'}
        Bf = [-diff(diff(N,x),x)'
             -diff(diff(N,y),y)'
             -2*diff(diff(N,x),y)']; % strain-displacement matrix
        Kf = double(int(int(Bf'*C*Bf,y,-dims.height/2,dims.height/2),x,-dims.width/2,dims.width/2));
        Ks = 0;
    % Mindlin
    case 'MB4'
        Bf = [diff(N(2,:),x)
              diff(N(3,:),x)
              diff(N(2,:),y) + diff(N(3,:),x)]; % flexural strain-displacement matrix
        Bs = [diff(N(1,:),x)
              diff(N(1,:),y)] + N(2:3,:); % shear strain-displacement matrix
        Kf = double(int(int(Bf'*C*Bf,y,-dims.height/2,dims.height/2),x,-dims.width/2,dims.width/2));
        Ks = double(int(int(Bs'*C*Bs,y,-dims.height/2,dims.height/2),x,-dims.width/2,dims.width/2));
end