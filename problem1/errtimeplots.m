%{

    from n = 3 to 1000, generate random equality constrained quadratic problem.

    Test the following solvers:
    - LDL dense
    - LU dense
    - Null-space
    - Range-space
    - LDL sparse
    - LU sparse

    Check their error and time taken to solve the problem.

    Make plots of error vs n, and time vs n for each solver.

%}

% Set up
fp = '/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/V1';
n = 3:10:1000;
num_trials = 10;
save_fig = false;
save_data = false;

% Define the colors
colors = [
    0, 0.4470, 0.7410;  % Blue
    0.8500, 0.3250, 0.0980;  % Red
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4940, 0.1840, 0.5560;  % Purple
    0.4660, 0.6740, 0.1880;  % Green
    0.3010, 0.7450, 0.9330;  % Cyan
];

% Initialize error and time arrays

error_LDLdense = zeros(length(n), num_trials);
error_LUdense = zeros(length(n), num_trials);
error_nullspace = zeros(length(n), num_trials);
error_rangespace = zeros(length(n), num_trials);
error_LDLsparse = zeros(length(n), num_trials);
error_LUsparse = zeros(length(n), num_trials);

time_LDLdense = zeros(length(n), num_trials);
time_LUdense = zeros(length(n), num_trials);
time_nullspace = zeros(length(n), num_trials);
time_rangespace = zeros(length(n), num_trials);
time_LDLsparse = zeros(length(n), num_trials);
time_LUsparse = zeros(length(n), num_trials);

% Run trials
sparsity = 15; % 15% sparsity

for i = 1:length(n)
    for j = 1:num_trials
        disp(['n = ', num2str(n(i)), ', trial = ', num2str(j)])
        % Generate random ECQP
        [H, g, A, b, x, lambda] = randECQP(n(i), sparsity);

    
        % Solve using dense LDL 
        tic;
        [x_sol, ~] = EqualityQPSolver(H, g, A, b, 'LDLdense');
        time_LDLdense(i, j) = toc;
        error_LDLdense(i, j) = norm(x - x_sol);
        
        % Solve using dense LU
        tic;
        [x_sol, ~] = EqualityQPSolver(H, g, A, b, 'LUdense');
        time_LUdense(i, j) = toc;
        error_LUdense(i, j) = norm(x - x_sol);
        
        % Solve using Null-space
        tic;
        [x_sol, ~] = EqualityQPSolver(H, g, A, b, 'nullspace');
        time_nullspace(i, j) = toc;
        error_nullspace(i, j) = norm(x - x_sol);
        
        % Solve using Range-space
        tic;
        [x_sol, ~] = EqualityQPSolver(H, g, A, b, 'rangespace');
        time_rangespace(i, j) = toc;
        error_rangespace(i, j) = norm(x - x_sol);
        
        % Solve using sparse LDL
        tic;
        [x_sol, ~] = EqualityQPSolver(H, g, A, b, 'LDLsparse');
        time_LDLsparse(i, j) = toc;
        error_LDLsparse(i, j) = norm(x - x_sol);
        
        % Solve using sparse LU
        tic;
        [x_sol, ~] = EqualityQPSolver(H, g, A, b, 'LUsparse');
        time_LUsparse(i, j) = toc;
        error_LUsparse(i, j) = norm(x - x_sol);
    end
end

% end

% Calculate mean error and time

mean_error_LDLdense = mean(error_LDLdense, 2);
mean_error_LUdense = mean(error_LUdense, 2);
mean_error_nullspace = mean(error_nullspace, 2);
mean_error_rangespace = mean(error_rangespace, 2);
mean_error_LDLsparse = mean(error_LDLsparse, 2);
mean_error_LUsparse = mean(error_LUsparse, 2);

mean_time_LDLdense = mean(time_LDLdense, 2);
mean_time_LUdense = mean(time_LUdense, 2);
mean_time_nullspace = mean(time_nullspace, 2);
mean_time_rangespace = mean(time_rangespace, 2);
mean_time_LDLsparse = mean(time_LDLsparse, 2);
mean_time_LUsparse = mean(time_LUsparse, 2);

% Plot error vs n 

figure;
hold on;
plot(n, mean_error_LDLdense, 'Color', colors(1, :));
plot(n, mean_error_LDLsparse, 'Color', colors(2, :));
plot(n, mean_error_LUdense, 'Color', colors(3, :));
plot(n, mean_error_LUsparse, 'Color', colors(4, :));
plot(n, mean_error_nullspace, 'Color', colors(5, :));
plot(n, mean_error_rangespace, 'Color', colors(6, :));

xlabel('n');
ylabel('$\Vert x - x^* \Vert$', 'Interpreter', 'latex');
legend('LDL', 'LDL sparse', 'LU', 'LU sparse', 'Null', 'Range');

if save_fig
    save_figure(gcf, 'errorplot2');
end
hold off;
% Plot time vs n 
figure;
hold on;
plot(n, mean_time_LDLdense,  'Color', colors(1, :));
plot(n, mean_time_LDLsparse, 'Color', colors(2, :));
plot(n, mean_time_LUdense, 'Color', colors(3, :));
plot(n, mean_time_LUsparse, 'Color', colors(4, :));
plot(n, mean_time_nullspace, 'Color', colors(5, :));
plot(n, mean_time_rangespace, 'Color', colors(6, :));
xlabel('n');
ylabel('Time (s)');
legend('LDL', 'LDL sparse', 'LU', 'LU sparse', 'Null', 'Range');
% saveas(gcf, '/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/Figures/timeplot.png');
if save_fig
    save_figure(gcf, 'timeplot');
end
% Save data (for future use)


% Errors
if save_data
    % save data to fp + 'error_LDLdense.mat'
    save(fullfile(fp, 'error_LDLdense.mat'), 'error_LDLdense');
    save(fullfile(fp, 'error_LDLsparse.mat'), 'error_LDLsparse');
    save(fullfile(fp, 'error_LUdense.mat'), 'error_LUdense');
    save(fullfile(fp, 'error_LUsparse.mat'), 'error_LUsparse');
    save(fullfile(fp, 'error_nullspace.mat'), 'error_nullspace');
    save(fullfile(fp, 'error_rangespace.mat'), 'error_rangespace');
    % Time
    save(fullfile(fp, 'time_LDLdense.mat'), 'time_LDLdense');
    save(fullfile(fp, 'time_LDLsparse.mat'), 'time_LDLsparse');
    save(fullfile(fp, 'time_LUdense.mat'), 'time_LUdense');
    save(fullfile(fp, 'time_LUsparse.mat'), 'time_LUsparse');
    save(fullfile(fp, 'time_nullspace.mat'), 'time_nullspace');
    save(fullfile(fp, 'time_rangespace.mat'), 'time_rangespace');

end



function save_figure(figure_handle, filename, resolution, width, height)
    % Save a MATLAB figure with specified attributes
    %
    % Parameters:
    %   figure_handle - Handle to the figure to be saved (default: current figure)
    %   filename - Name of the file to save (default: 'figure')
    %   resolution - Resolution in dpi (default: 300)
    %   width - Width of the figure in inches (default: 6 inches)
    %   height - Height of the figure in inches (default: 4 inches)
    %
    % Example usage:
    %   save_figure(gcf, 'my_figure', 600, 8, 6);

    if nargin < 1 || isempty(figure_handle)
        figure_handle = gcf;
    end
    if nargin < 2 || isempty(filename)
        filename = 'figure';
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
    print(figure_handle, filename, '-dpng', ['-r', num2str(resolution)]);
    fprintf('Figure saved as %s.png at %d dpi\n', filename, filename, filename, resolution);
end
