function plotConvergence(iter, values, obj)
% Plot convergence plots.
% 'iter' is a vector containing the number of iterations.
% 'values' is a vector containing the values to plot.
% 'obj' is a caracter representing the type of convergence plot.

switch obj
    case 'c'
        % plot convergence in terms of compliance
        figure(1);
        plot(iter, values);
        title('Convergence of the objective function: Compliance');
        xlabel('Iterations'); ylabel('Compliance');
        drawnow;
    case 'f'
        % plot convergence in terms of the lowest eigenfrequency
        figure(1);
        plot(iter, values);
        names = '';
        for i=1:size(values, 1)
            names = strcat(names, sprintf('%.f°,', i));
        end
        legend(strsplit(names, ','));
        title(sprintf('Convergence of the objective function: %.f° frequency', size(values, 1)-3));
        xlabel('Iterations'); ylabel(sprintf('%.f° frequency', size(values, 1)));
        drawnow;
    case 'x'
        % plot convergence in terms of density change
        figure(2);
        plot(iter, values);
        title('Convergence of the design variable: Density');
        xlabel('Iterations'); ylabel('Density change');
        drawnow;
end
end