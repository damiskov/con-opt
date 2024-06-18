%{
    Driver code for problem 4 - SQP for Himmelblau's test problem
%}

step_comp = 'trust_region'; % 'line_search' or 'trust_region
hessian_method = 'hessian'; % 'BFGS' or 'hessian';
% Initial guesses
x1 =  [-2; 0.8];
x2 = [-5; 0];
x3 =  [1; 0];
x4 = [1; 5];



start_points_eqcon = {[0; 0],[5; 5], [-2; 0], [-5; -2], [5; 0]};
% Lower and upper bounds
x_l = [-5; -5];
x_u = [5; 5]; 
% Lower and upper bounds for inequality constraints
g_l = 0;
g_u =  1e6; % Large number for upper bound (inf)
objective = @himmelblau;
% Inequality constraint
g = @himmelblauIneq;

equality_constraint = false; % Set to true to include equality constraint

if equality_constraint
    start_points = start_points_eqcon;
    fname = 'Himmelblau_eqcon_';
    h = @himmelblauEq;
else
    start_points = {x1, x2, x3, x4};
    fname = 'Himmelblau_no_eqcon_';
    h = @noEqCon;
end


for i = 1:length(start_points)
    x0 = start_points{i};
    disp(x0)

    if strcmp(step_comp, 'line_search')
        if strcmp(hessian_method, 'BFGS')
            [x, f, exitflag, iterations, lambda, DxL, H, KKT] = SQPlinesearch(objective, x0, x_l, x_u, g_l, g_u, h, g);
            plotHBIter(x', equality_constraint);
            save_figure(gcf, ['figures/',fname,'_BFGS_', step_comp, '_start_', num2str(i)], 300, 6, 4);
            save(['data/',fname,'_BFGS_', step_comp, '_start_', num2str(i), '.mat'], 'x', 'f', 'DxL', 'KKT');
        else
            [x, f, exitflag, iterations, lambda, DxL, H, KKT] = SQPlinesearch(objective, x0, x_l, x_u, g_l, g_u, h, g, @himmelblauHessian);
            plotHBIter(x', equality_constraint);
            save_figure(gcf, ['figures/',fname,'_hessian_', step_comp, '_start_', num2str(i)], 300, 6, 4);

            save(['data/',fname,'_hessian_', step_comp, '_start_', num2str(i), '.mat'], 'x', 'f', 'DxL', 'KKT');
        end
    end


    if strcmp(step_comp, 'trust_region')
        if strcmp(hessian_method, 'BFGS')
            [x, f, exitflag, iterations, lambda, DxL, H, KKT] = SQPtrustregion(objective, x0, x_l, x_u, g_l, g_u, h, g);

            plotHBIter(x', equality_constraint);
            save_figure(gcf, ['figures/',fname,'_BFGS_', step_comp, '_start_', num2str(i)], 300, 6, 4);

            save(['data/',fname,'_BFGS_', step_comp, '_start_', num2str(i), '.mat'], 'x', 'f', 'DxL', 'KKT');
        else
            [x, f, exitflag, iterations, lambda, DxL, H, KKT] = SQPtrustregion(objective, x0, x_l, x_u, g_l, g_u, h, g, @himmelblauHessian);

            plotHBIter(x', equality_constraint);
            save_figure(gcf, ['figures/',fname,'hessian_', step_comp, '_start_', num2str(i)], 300, 6, 4);

            save(['data/',fname,'_hessian_', step_comp, '_start_', num2str(i), '.mat'], 'x', 'f', 'DxL', 'KKT');
        end
    end

end




function [f, df, d2f] = himmelblau(x)
    f = (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
    if nargout > 1
        df = [4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7); 2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7)];
        if nargout > 2
            d2f = [12*x(1)^2 + 4*x(2) - 42, 4*(x(1) + x(2)); 4*(x(1) + x(2)), 4*x(1) + 12*x(2)^2 - 26];
        end
    end
end


function [h, dh, d2h] =  himmelblauEq(x)
    %{
        Equality constraint for Himmelblau's test problem.

        y = 2/3 x
    %}

    h = [(2/3)*x(1) - x(2)];
    dh = [2/3;-1];
    d2h = [[0 0;0 0], [0 0;0 0]];
end

function [h, dh, d2h] = noEqCon(x)
    h = 0;
    dh = [0;0];
    d2h = [[0 0;0 0], [0 0;0 0]];
end


function [g, dg, d2g] = himmelblauIneq(x)
    g = [(x(1)+2)^2 - x(2); -4*x(1) + 10*x(2)];
    dg = [2*(x(1)+2), -4; -1, 10];
    d2g = [[2 0;0 0], [0 0;0 0], [-2 0;0 0], [0 0;0 0]];
end

function H = himmelblauHessian(x)
    %{
        Hessian of the objective function for our Himmelblau test problem.
    %}

    H = [12*x(1)^2 + 4*x(2) - 42, 4*x(1) + 4*x(2);
     4*x(1) + 4*x(2), 4*x(1) + 12*x(2)^2 - 26];
end



function plotHBIter(x_iter, eqcon)

    % Plotting iterations on himmelblau contour plot along with feasible region and equality constraints
    % x: data from iterations


    x1 = linspace(-6, 6, 100);
    x2 = linspace(-6, 6, 100);
    [X1, X2] = meshgrid(x1, x2);


    Z = zeros(size(X1));
    for i = 1:size(X1, 1)
        for j = 1:size(X1, 2)
            Z(i, j) = himmelblau([X1(i, j); X2(i, j)]);
        end
    end

    % Plot the contour
    figure;
    v = [-10:0.5:0, 0:1:10 10:5:100 100:20:200, 200:50:400];
    contour(X1, X2, Z, v);

    % Shading the feasible region

    yc1 = (x1+2).^2;
    yc2 = 0.4.*x1;

    % Adding equality constraint: y = x, to plot

    % Adding constraints+feasible regions to the plot

    hold on

        % Equality constraint that goes through the feasible minimum, but not through the origin

        if eqcon 
            yeq = 2/3*x1;

            plot(x1, yeq, 'LineWidth', 1, 'Color', 'black');
        end
    
        
        % Inequality constraints
        fill(x1,yc1,[0.7 0.7 0.7],'facealpha',0.5)
        fill([x1 x1(end) x1(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.5)

        % Box Constraints
        % Shading area below lower bound
        fill([x1(1), x1(1), x1(end), x1(end)], [-6, -5, -5, -6], [0.7 0.7 0.7], 'facealpha', 0.5)
        % Shading area above upper bound
        fill([x1(1), x1(1), x1(end), x1(end)], [5, 6, 6, 5], [0.7 0.7 0.7], 'facealpha', 0.5)
        % Shading area to the left of the left bound
        fill([-6, -5, -5, -6], [x2(1), x2(1), x2(end), x2(end)], [0.7 0.7 0.7], 'facealpha', 0.5)
        % Shading area to the right of the right bound
        fill([5, 6, 6, 5], [x2(1), x2(1), x2(end), x2(end)], [0.7 0.7 0.7], 'facealpha', 0.5)
        
    hold off

    xlim([-6 6])
    ylim([-6 6])
    colorbar
    xlabel('x_1');
    ylabel('x_2');

    hold on

    x = x_iter(1, :);
    y = x_iter(2, :);

    
    plot(x(1), y(1), '*', 'Color', [0.4660, 0.6740, 0.1880], 'MarkerSize', 10, 'LineWidth', 2)


    plot(x(1:end), y(1:end), '-x', 'Color', [0.8500, 0.3250, 0.0980], 'MarkerSize', 10, 'LineWidth', 2)

    % Mark final point
    % Blue RBG: [0, 0.4470, 0.7410]
    plot(x(end), y(end), '+', 'Color', [0, 0.4470, 0.7410], 'MarkerSize', 10, 'LineWidth', 2)

end


function save_figure(figure_handle, filepath, resolution, width, height)
    % Save a MATLAB figure with specified attributes
    %
    % Parameters:
    %   figure_handle - Handle to the figure to be saved (default: current figure)
    %   filepath - path to save file (default: 'figure')
    %   resolution - Resolution in dpi (default: 300)
    %   width - Width of the figure in inches (default: 6 inches)
    %   height - Height of the figure in inches (default: 4 inches)
    %
    % Example usage:
    %   save_figure(gcf, 'my_figure', 600, 8, 6);

    if nargin < 1 || isempty(figure_handle)
        figure_handle = gcf;
    end
    if nargin < 2 || isempty(filepath)
        filepath = 'figure';
    end
    if nargin < 3 || isempty(resolution)
        resolution = 300;
    end
    if nargin < 4 || isempty(width)
        width = 6;
    end
    if nargin < 5 || isempty(height)
        height = 4;
    end

    % Set the paper size
    set(figure_handle, 'PaperUnits', 'inches');
    set(figure_handle, 'PaperPosition', [0 0 width height]);
    set(figure_handle, 'PaperSize', [width height]);
    set(figure_handle, 'Renderer', 'painters');
    print(figure_handle, filepath, '-dpng', ['-r', num2str(resolution)]);
    fprintf('Figure saved as %s.png at %d dpi\n', filepath, resolution);
end

