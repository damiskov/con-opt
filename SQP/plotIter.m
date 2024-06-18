function plotIter(x_iter, eqcon)

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