% Script to read CSV file and plot multiple pressure transducer readings
% The first column (Timestamp) will be the x-axis
% Usage: Simply run this script and select your CSV file when prompted

% Prompt user to select a CSV file
[fileName, filePath] = uigetfile('*.csv', 'Select the CSV file');
if fileName == 0
    disp('No file selected. Exiting script.');
    return;
end

% Read the CSV file
fullPath = fullfile(filePath, fileName);
data = readtable(fullPath);

% Find PT columns (columns with "Pressure" in the name)
ptColumns = find(contains(columnNames, 'Pressure'));
disp(['Found ', num2str(length(ptColumns)), ' pressure transducer columns:']);
for i = 1:length(ptColumns)
    disp(['  ', num2str(i), '. ', columnNames{ptColumns(i)}]);
end

numPTsToPlot = length(ptColumns);

% Get the timestamp column (first column)
timestampCol = columnNames{1};
timestamps = data.(timestampCol);

% Create a figure
figure;
hold on;

% Color map for multiple lines
colors = lines(numPTsToPlot);

% Plot each selected PT column
legendEntries = cell(numPTsToPlot, 1);
for i = 1:numPTsToPlot
    colIdx = ptColumns(i);
    colName = columnNames{colIdx};
    plot(timestamps, data.(colName), 'LineWidth', 1.5, 'Color', colors(i,:));
    legendEntries{i} = colName;
end

% Add labels and legend
xlabel('Timestamp', 'Interpreter', 'none');
ylabel('Pressure', 'Interpreter', 'none');
title('Pressure Transducer Readings vs Time', 'Interpreter', 'none');
legend(legendEntries, 'Location', 'best', 'Interpreter', 'none');
grid on;

% % Adjust figure appearance
% set(gcf, 'Position', [100, 100, 1200, 800]);
% set(gca, 'FontSize', 12);

% If timestamps are datetime objects, rotate x-axis labels for better readability
if isdatetime(timestamps)
    xtickangle(45);
end

hold off;
disp('Plot created successfully!');