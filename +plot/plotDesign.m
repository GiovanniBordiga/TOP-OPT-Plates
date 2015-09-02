function plotDesign(x)
% Plot the density field.
% 'x' is a matrix representing the density field on the plate.

figure(3);                      % focus the corresponding figure
colormap(flipud(gray));         % set a gray-scale colormap
imagesc(x);                     % display 'x'
grid;                           % show the grid
colorbar;                       % display colorbar
axis equal; axis tight;         % set axis to match plate dimensions
title('Density field');         % set title
xlabel('x'); ylabel('y');       % set axis labels
drawnow;
end