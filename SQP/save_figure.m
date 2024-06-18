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