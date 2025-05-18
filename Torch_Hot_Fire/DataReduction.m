% Use file selection dialog to choose the CSV file
[filename, filepath] = uigetfile({'*.csv', 'CSV Files (*.csv)'; '*.*', 'All Files (*.*)'}, 'Select CSV Data File');

% Check if user canceled the dialog
if isequal(filename, 0) || isequal(filepath, 0)
    disp('File selection canceled');
    return;
end

% Construct full file path
fullFilePath = fullfile(filepath, filename);

try
    % Read the CSV file with preserved variable names
    opts = detectImportOptions(fullFilePath);
    opts.VariableNamingRule = 'preserve'; % Keep original column names
    data = readtable(fullFilePath, opts);
    
    % Display variable names and first few rows for debugging
    disp('Available column names:');
    disp(data.Properties.VariableNames);
    disp('First few rows of data:');
    disp(data(1:3,:));
    
    % Convert timestamps and create proper time vector
    baseTime = datetime(data.Timestamp(1));
    numPoints = height(data);
    
    % Create a proper time vector with 10 points per second
    timeVector = baseTime + seconds((0:(numPoints-1))/10);
    
    % Create a figure
    figure('Position', [100, 100, 1000, 600]);
    
    % Plot pressure values
    plot(timeVector, data{:, 'PT-TI-01 Pressure'}, 'b-', 'LineWidth', 2);
    hold on;
    plot(timeVector, data{:, 'PT-O2-05 Pressure'}, 'g-', 'LineWidth', 2);
    plot(timeVector, data{:, 'PT-H2-03 Pressure'}, 'r-', 'LineWidth', 2);
    
    % Using a secondary Y-axis for N2 pressure which has much higher values
    yyaxis right;
    plot(timeVector, data{:, 'PT-N2-04 Pressure'}, 'm-', 'LineWidth', 2);
    ylabel('N2 Pressure');
    yyaxis left;
    
    grid on;
    % title('Pressure Measurements with Solenoid State Changes');
    xlabel('Time');
    ylabel('Pressure');
    
    % Get the current axis limits for drawing vertical lines
    yyaxis left;
    ylim1 = ylim;
    
    % Find state changes and add vertical lines
    stateVars = {'SN-H2-01 State', 'SN-O2-01 State', 'SN-N2-01 State', 'Spark Plug State'};
    stateColors = {'r', 'g', 'm', 'k'};
    stateLabels = {'H2', 'O2', 'N2', 'Spark'};
    
    % Process each state variable
    for i = 1:length(stateVars)
        % Get state data and convert to numeric if it's text/cell
        stateData = data{:, stateVars{i}};
        
        % Check data type and convert if needed
        if iscell(stateData)
            % If it's a cell array (text), convert 'True'/'False' to numeric
            numericState = zeros(size(stateData));
            for k = 1:length(stateData)
                if strcmpi(stateData{k}, 'True') || strcmpi(stateData{k}, 'true') || ...
                   strcmpi(stateData{k}, 'T') || strcmpi(stateData{k}, '1')
                    numericState(k) = 1;
                else
                    numericState(k) = 0;
                end
            end
            stateData = numericState;
        elseif islogical(stateData)
            % If it's already logical, convert to double
            stateData = double(stateData);
        end
        
        % Find indices where state changes
        stateChanges = find(diff(stateData) ~= 0) + 1;
        
        % Only proceed if there are state changes
        if ~isempty(stateChanges)
            % Get state at each change point
            stateValues = stateData(stateChanges);
            
            % For legend entries
            stateOn = line(NaN, NaN, 'Color', stateColors{i}, 'LineStyle', '-', 'LineWidth', 2);
            
            % Plot markers for each state change
            for j = 1:length(stateChanges)
                idx = stateChanges(j);
                currentState = stateValues(j);
                
                % Draw a vertical line at each state change
                line([timeVector(idx), timeVector(idx)], ylim1, 'Color', stateColors{i}, 'LineStyle', '-', 'LineWidth', 1.5);
                
                % Add a text label above the plot to indicate which solenoid changed
                if currentState > 0
                    % State turning ON
                    text(timeVector(idx), ylim1(2)*1.05, [stateLabels{i} ' ON'], 'Color', stateColors{i}, ...
                         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                else
                    % State turning OFF
                    text(timeVector(idx), ylim1(2)*1.05, [stateLabels{i} ' OFF'], 'Color', stateColors{i}, ...
                         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                end
            end
        end
    end
    
    % Create legend with both pressure and state markers
    pressureLegend = {'PT-TI-01', 'PT-O2-05', 'PT-H2-03', 'PT-N2-04'};
    stateLegend = {'H2 State Change', 'O2 State Change', 'N2 State Change', 'Spark State Change'};
    
    % Create dummy line objects for the legend
    hLines = zeros(length(stateVars), 1);
    for i = 1:length(stateVars)
        hLines(i) = line(NaN, NaN, 'Color', stateColors{i}, 'LineStyle', '-', 'LineWidth', 1.5);
    end
    
    % Combine legends
    legend([pressureLegend, stateLegend], 'Location', 'best');
    
    % Set format for time display on x-axis to show milliseconds
    datetick('x', 'HH:MM:SS.FFF', 'keepticks');
    
    % % Set the title with the filename
    % sgtitle(['Pressure Monitoring with Solenoid State Changes - ' filename]);
    
catch ME
    % Display any errors
    errordlg(['Error processing file: ' ME.message], 'Error');
    disp(['Error details: ' ME.getReport]);
end