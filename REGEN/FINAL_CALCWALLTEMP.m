function [flowTemp,flowVelocity,flowPressure, T_wgFinal, finEfficiency, Qdot, finQdot, T_wl_Array, h_l_Array, h_g_Array, vonMises, sigma_long, sigma_circ, sigma_rad] = FINAL_CALCWALLTEMP(converge_index, throat_index, channelHeight, channelWidth, wallThickness, chamberDiameterArray, heightStepNumber, newFluidProperties, inputValues)
    %% Channel Dimensions
   
    channelWidths = [linspace(channelWidth(1),channelWidth(2),throat_index),linspace(channelWidth(2),channelWidth(3),converge_index-throat_index),(channelWidth(3)*ones(1,(heightStepNumber - converge_index)))];
    channelHeights = [linspace(channelHeight(1),channelHeight(2),throat_index),linspace(channelHeight(2),channelHeight(3),converge_index-throat_index),(channelHeight(3)*ones(1,(heightStepNumber - converge_index)))];
    wallThicknesses = [linspace(wallThickness(1),wallThickness(2),throat_index),linspace(wallThickness(2),wallThickness(3),converge_index-throat_index),(wallThickness(3)*ones(1,(heightStepNumber - converge_index)))];

    
    %% Call Heat Transfer Solver
    [T_wgFinal, flowTemp, flowPressure, flowVelocity, finEfficiency, Qdot, finQdot, T_wl_Array, h_l_Array, h_g_Array] = FINAL_SOLVER(chamberDiameterArray, heightStepNumber, newFluidProperties, channelWidths, channelHeights, wallThicknesses,inputValues);
    
    %% Structural Analysis
    chamberPressure = newFluidProperties(:,9)*100000; % Chamber Pressure (Pa)
    deltaP = flowPressure - chamberPressure';
    [vonMises,sigma_long, sigma_circ, sigma_rad, structuresOutput] = FINAL_STRUCTURES(T_wl_Array, T_wgFinal, wallThicknesses, channelWidths, chamberPressure', chamberDiameterArray', deltaP, Qdot, inputValues);
    

    %% Plotter: Structures
    
    % Circ stress
        axialDist = newFluidProperties(:,1)* 39.3701;
        
    
        figure;
        hold on; box on;
        
        % Plot with contrasting, muted colors (blue, orange, gray)
        p1 = plot(axialDist, sigma_circ,  '-', 'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]); % blue
        p2 = plot(axialDist, structuresOutput(2,:), '-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]); % orange
        p3 = plot(axialDist, structuresOutput(1,:),  '-', 'LineWidth', 1.8, 'Color', [0.50 0.50 0.50]); % gray
        
        % Axis labels
        xlabel('Location [in]', 'FontSize', 12);
        ylabel('Circumferential stress [MPa]', 'FontSize', 12);
        
        % Grid and ticks
        grid on;
        grid minor;
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
                 'FontSize', 11, 'LineWidth', 1);
        xlim([0, axialDist(1)]);
        
        ylim([0, max(sigma_circ)*1.1]);
        % Legend
        legend([p1 p2 p3], {'Total Stress', 'Pressure Load', 'Thermal Load'}, ...
               'Location', 'best', 'Box', 'off');


     % Channel Dimensions
        figure;
        hold on; box on;
        
        % Plot with contrasting, muted colors (blue, orange, gray)
        h1 = plot(axialDist, channelHeights, '-', 'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]); % blue
        h2 = plot(axialDist, channelWidths,  '-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]); % orange
        h3 = plot(axialDist, wallThicknesses, '-', 'LineWidth', 1.8, 'Color', [0.50 0.50 0.50]); % gray
        
        % Axis labels
        xlabel('Location [in]', 'FontSize', 12);
        ylabel('Channel Dimensions [mm]', 'FontSize', 12);
        
        % Grid and ticks
        grid on;
        grid minor;
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
                 'FontSize', 11, 'LineWidth', 1);
        xlim([0, axialDist(1)]);
        
        ylim([0, max([channelWidths,channelHeights,wallThicknesses])*1.1]);
        % Legend
        legend([h1 h2 h3], {'Channel Height', 'Channel Width', 'Wall Thickness'}, ...
               'Location', 'best', 'Box', 'off');
    
    % Von Mises vs Yield
        %6061ram2:
        %yieldStrength = (-(1.57*10^-5)*(T_wgFinal.^3) + (0.0139 * (T_wgFinal.^2)) - 4.21*(T_wgFinal) + 762)*1000000;
        %UTS =  (-(1.31*10^-5)*(T_wgFinal.^3) + (0.0113 * (T_wgFinal.^2)) - 3.45*(T_wgFinal) + 719)*1000000;
        %7075:
        yieldStrength1 = (-19756 + (188.4887580831*T_wgFinal) -(.63101218051*(T_wgFinal.^2)) + (.0008997028*(T_wgFinal.^3)) -(0.0000004663*(T_wgFinal.^4)))*1000000;
        UTS1 = (-20399.3454015107 + (196.2160186167*T_wgFinal)  -(0.6601708400*(T_wgFinal.^2)) + (.0009446715*(T_wgFinal.^3)) -(.0000004910*(T_wgFinal.^4)))*1000000;
        
        yieldStrength2 = (-19756 + (188.4887580831*T_wl_Array) -(.63101218051*(T_wl_Array.^2)) + (.0008997028*(T_wl_Array.^3)) -(0.0000004663*(T_wl_Array.^4)))*1000000;
        UTS2 = (-20399.3454015107 + (196.2160186167*T_wl_Array)  -(0.6601708400*(T_wl_Array.^2)) + (.0009446715*(T_wl_Array.^3)) -(.0000004910*(T_wl_Array.^4)))*1000000;
         

        % yieldStrength1 = 297.4 - 0.312*T_wgFinal + 3.175*(10^-4).*(T_wgFinal.^2) - 2.506*(10^-7).*(T_wgFinal.^3); % GRCOP42
        % UTS1 = 664.5 - 0.987*T_wgFinal + 3.58*(10^-4).*(T_wgFinal.^2);
        % 
        % yieldStrength2 = 297.4 - 0.312*T_wl_Array + 3.175*(10^-4).*(T_wl_Array.^2) - 2.506*(10^-7).*(T_wl_Array.^3); % GRCOP42
        % UTS2 = 664.5 - 0.987*T_wl_Array + 3.58*(10^-4).*(T_wl_Array.^2);
        % 
        yieldStrength = (yieldStrength1+yieldStrength2)/2;
        UTS = (UTS1+UTS2)/2;
        % Create figure
        figure;
        box on;
        
        % LEFT Y-AXIS: Strength & Stress [MPa]
        yyaxis left
        hold on;
        
        % Muted, contrasting colors: blue, orange, gray
        h1 = plot(axialDist * 39.3701, vonMises/1e6,'-', 'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]); % blue
        h2 = plot(axialDist * 39.3701, yieldStrength/1e6,'-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]); % orange
        h3 = plot(axialDist * 39.3701, UTS/1e6, '-', 'LineWidth', 1.8, 'Color', [0.50 0.50 0.50]); % gray
        
        ylabel('Strength & Stress [MPa]', 'FontSize', 12);
        
        
        allY_left = [UTS,vonMises]/1e6;   % divide by 1e6 here too
        ylim([0, max(allY_left)*1.1]);

        
        
        %% RIGHT Y-AXIS: Contour (Radial Distance [m])
        yyaxis right
        
        h4 = plot(axialDist * 39.3701, chamberDiameterArray/2 * 39.3701, 'k--', 'LineWidth', 1.8);
        ylabel('Radial Distance [in]', 'FontSize', 12);
        
        % --- Control RIGHT Y limits: 0 to next grid mark above max Contour ---
        
      
        ylim([0, max(chamberDiameterArray)*1.1 * 39.3701]);
       
        
        %% Common X-axis
        xlabel('Axial Distance [in]', 'FontSize', 12);
        
        % X limits: start and end at data bounds
        xlim([0,axialDist(1) * 39.3701]);
        
        %% Grid, ticks, formatting
        grid on;
        grid minor;
        
        ax = gca;
        set(ax, 'XMinorTick', 'on', ...
                'YMinorTick', 'on', ...
                'FontSize', 11, ...
                'LineWidth', 1);
        
        % Make both y-axes use dark text
        ax.YColor = [0 0 0];
        
        %% Legend
        legend([h1 h2 h3 h4], ...
               {'Effective Stress', 'Yield Strength', 'Ultimate Tensile Strength', 'Contour'}, ...
               'Location', 'best', 'Box', 'off');    

        


    % Temperatures Plot
        
    figure;

        % Plot with muted, contrasting colors (blue, orange, gray)
        plot(axialDist, flowTemp, '-', 'LineWidth', 1.8, 'Color', [0 0.45 0.74]);  % blue
        hold on;
        plot(axialDist, T_wgFinal, '-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]);   % orange
        plot(axialDist, T_wl_Array, '-', 'LineWidth', 1.8, 'Color', [0.4 0.4 0.4]);     % gray
       
        hold off;
        
        % Axes labels
        xlabel('Location [in]', 'FontSize', 12);
        ylabel('Temperature [K]', 'FontSize', 12);
        ylim([0,max([T_wl_Array,T_wgFinal,flowTemp])*1.1])
        % Title
        title('Fluid and Wall Temperatures along Chamber', ...
              'FontSize', 14, 'FontWeight', 'bold');
        
        % Legend
        legend({'Flow Temp', 'Gas Hotwall Temp', 'Liquid Hotwall Temp'}, ...
               'Location', 'best', 'Box', 'off');
        
        % Grid and ticks formatting
        grid on;
        grid minor;
        set(gca, 'LineWidth', 1, 'FontSize', 11, 'TickDir', 'out', 'Box', 'on');

        yyaxis right
        
        plot(axialDist, chamberDiameterArray/2 * 39.3701, 'k--', 'LineWidth', 1.8);
        ylabel('Radial Distance [in]', 'FontSize', 12);
        
        % --- Control RIGHT Y limits: 0 to next grid mark above max Contour ---
        
      
        ylim([0, max(chamberDiameterArray)*1.1 * 39.3701]);
        legend({'Flow Temp', 'Gas Hotwall Temp', 'Liquid Hotwall Temp' ,'Contour'}, ...
               'Location', 'best', 'Box', 'off');

    %Flow Pressure and Velocity
        figure;

        % Plot with muted, contrasting colors (blue, orange, gray)
        plot(axialDist, flowPressure/6895, '-', 'LineWidth', 1.8, 'Color', [0 0.45 0.74]);  % blue
        hold on;
        plot(axialDist, chamberPressure/6895, '-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]);   % orange
        
        hold off;
        
        % Axes labels
        xlabel('Location [in]', 'FontSize', 12);
        ylabel('Pressure [psi]', 'FontSize', 12);
        ylim([0,max([flowPressure,chamberPressure'])*1.1/6895])
        % Title
        title('Coolant & Exhaust Pressures & Velocity', ...
              'FontSize', 14, 'FontWeight', 'bold');
        
        % Legend
       
        
        % Grid and ticks formatting
        grid on;
        grid minor;
        set(gca, 'LineWidth', 1, 'FontSize', 11, 'TickDir', 'out', 'Box', 'on');

        yyaxis right
        
        plot(axialDist, flowVelocity*3.28084,'-', 'LineWidth', 1.8, 'Color', [0.4 0.4 0.4]);
        ylabel('Velocity (ft/s)', 'FontSize', 12);
        
        % --- Control RIGHT Y limits: 0 to next grid mark above max Contour ---
        
      
        ylim([0, max(flowVelocity)*1.1 * 3.28084]);
        legend({'Flow Pressure', 'Chamber Pressure','Flow Velocity'}, ...
               'Location', 'best', 'Box', 'off');
        

    
    % h_l, h_g
        
        figure;

        % Plot with muted, contrasting colors (blue, orange, gray)
        plot(axialDist, h_g_Array, '-', 'LineWidth', 1.8, 'Color', [0 0.45 0.74]);  % blue
        hold on;
        plot(axialDist, h_l_Array, '-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]);   % orange
        
        hold off;
        
        % Axes labels
        xlabel('Location [in]', 'FontSize', 12);
        ylabel('[W/(M^2*K)]', 'FontSize', 12);
        ylim([0,max([h_g_Array,h_l_Array])*1.1])
        % Title
        title('Convective Heat Transfer Coefficients (CHTC)', ...
              'FontSize', 14, 'FontWeight', 'bold');
        
        
        
        % Grid and ticks formatting
        grid on;
        grid minor;
        set(gca, 'LineWidth', 1, 'FontSize', 11, 'TickDir', 'out', 'Box', 'on');

        yyaxis right
        
        plot(axialDist, chamberDiameterArray/2 * 39.3701, 'k--', 'LineWidth', 1.8);
        ylabel('Radial Distance (in)', 'FontSize', 12);
        
        % --- Control RIGHT Y limits: 0 to next grid mark above max Contour ---
        
      
        ylim([0,max(chamberDiameterArray)*1.1 * 39.3701]);
        % Legend
        legend({'Gas-Side CHTC', 'Liquidx-Side CHTC','Contour'}, ...
               'Location', 'best', 'Box', 'off');
    % fin_width, height, efficiency
        
end