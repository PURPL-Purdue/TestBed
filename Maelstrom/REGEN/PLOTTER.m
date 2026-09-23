%% Plotter for varioous things
% Set equations and arrays
temp = [298:805];
CTE = (0.007111*temp) + 21.48;
% Channel Dimensions
        figure;
        hold on; box on;
        
        % Plot with contrasting, muted colors (blue, orange, gray)
        h1 = plot(temp, CTE, '-', 'LineWidth', 4.8, 'Color', [0.00 0.45 0.74]); % blue
        %h2 = plot(temp, UTS1,  '-', 'LineWidth', 4.8, 'Color', [0.85 0.33 0.10]); % orange
        % h3 = plot(axialDist, wallThicknesses, '-', 'LineWidth', 1.8, 'Color', [0.50 0.50 0.50]); % gray
        
        % Axis labels
        xlabel('Temperature [K]', 'FontSize', 16);
        ylabel('CTE [\mum/m*\circK]', 'FontSize', 16);
        
        % Grid and ticks
        grid on;
        grid minor;
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
                 'FontSize', 16, 'LineWidth', 1);
        xlim([temp(1), temp(end)]);
        
        ylim([15, max(CTE)*1.1]);
        % Legend
        legend( { 'CTE'}, ...
               'Location', 'best', 'Box', 'off');