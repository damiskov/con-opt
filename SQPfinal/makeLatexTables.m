%{
    Goes through all mat files in data folder, converts them to tables and then to latex tables.
%}

makeTables('fmincon', 'fmincon');


function makeTables(datafolder, savefolder)
    % Get all mat files in data folder
    files = dir([datafolder,'/*.mat']);
    for file = files'
        fprintf('Converting %s to latex table\n', file.name);
        % Convert to table
        x = load([datafolder,'/', file.name]);
        % Convert to latex table
        x = struct2table(x);
        % Need to reverse columns
        x = x(:, end:-1:1);
        % Write to file
        table2latex(x, [savefolder,'/', file.name(1:end-4), '.tex']);
    end
end