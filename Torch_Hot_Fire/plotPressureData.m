function plotPressureData(csvFilePath)
    % Function to plot pressure data from a CSV file
    % Usage: plotPressureData('your_pressure_data.csv')
    
    % Read the CSV file
    data = readtable(csvFilePath);
    
    % Extract the timestamp column and convert to datetime
    timestamps = datetime(data.Timestamp);
    
    % Extract the pressure columns we want to plot (excluding PT-N2-04)
    pressureColumns = {'PT_TI_01 Pressure', 'PT_O2_05 Pressure', 'PT_H2_03 Pressure'};
    
    % Create column names with underscores instead of hyphens for MATLAB compatibility
    columnNames = strrep(data.Properties.VariableNames, '-', '_');
    data.Properties.VariableNames = columnNames;
    
    % Create the figure
    figure('Name', 'Pressure Sensor Data', 'Position', [100, 100, 900, 500]);
    
    % Plot the data
    hold on;
    
    % Plot PT-TI-01 Pressure in blue
    plot(timestamps, data.PT_TI_01_Pressure, 'b-', 'LineWidth', 2, 'DisplayName', 'PT-TI-01 Pressure');
    
    % Plot PT-O2-05 Pressure in red
    plot(timestamps, data.PT_O2_05_Pressure, 'r-', 'LineWidth', 2, 'DisplayName', 'PT-O2-05 Pressure');
    
    % Plot PT-H2-03 Pressure in green
    plot(timestamps, data.PT_H2_03_Pressure, 'g-', 'LineWidth', 2, 'DisplayName', 'PT-H2-03 Pressure');
    
    % Add title and labels
    title('Pressure Sensors Over Time', 'FontSize', 14);
    xlabel('Time', 'FontSize', 12);
    ylabel('Pressure', 'FontSize', 12);
    
    % Add legend
    legend('Location', 'best');
    
    % Add grid
    grid on;
    
    % Format plot
    ax = gca;
    ax.FontSize = 11;
    
    % Rotate x-axis labels if needed for better readability
    xtickangle(45);
    
    hold off;
end