function dCnew = filterSensitivity(nelx, nely, R, x, dC)
% Implement a mesh-indipendent sensitivity filter, using a linear weight function.
% 'nelx' and 'nely' are the number of element along the two dimensions.
% 'x' is a nely-by-nelx matrix representing the density field on the plate.
% 'R' is the filter's radius

% WARNING: this code is written assuming unit element's area.

dCnew = zeros(nely, nelx);
for elx = 1:nelx
    for ely = 1:nely
        sumw = 0;
        for i = max(elx-R, 1):min(elx+R, nelx)
            for j = max(ely-R, 1):min(ely+R, nely)
                w = max(0, R - sqrt((i-elx)^2 + (j-ely)^2));
                dCnew(ely, elx) = dCnew (ely, elx) + w*x(j, i)*dC(j, i);
                sumw = sumw + w;
            end
        end
        dCnew(ely, elx) = dCnew(ely, elx)/(x(ely, elx)*sumw);
    end
end
end

