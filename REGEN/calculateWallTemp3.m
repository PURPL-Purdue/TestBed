function [flowTemp,flowVelocity,flowPressure, T_wgFinal, finEfficiency, Qdot, finQdot, T_wl_Array, h_l_Array, h_g_Array, vonMises, sigma_long, sigma_circ, sigma_rad] = calculateWallTemp3(converge_index, throat_index, channelHeight, channelWidth, wallThickness, chamberDiameterArray, heightStepNumber, newFluidProperties, inputValues)
    %% Channel Dimensions
   
    channelWidths = [linspace(channelWidth(1),channelWidth(2),throat_index),linspace(channelWidth(2),channelWidth(3),converge_index-throat_index),(channelWidth(3)*ones(1,(heightStepNumber - converge_index)))];
    channelHeights = [linspace(channelHeight(1),channelHeight(2),throat_index),linspace(channelHeight(2),channelHeight(3),converge_index-throat_index),(channelHeight(3)*ones(1,(heightStepNumber - converge_index)))];
    wallThicknesses = [linspace(wallThickness(1),wallThickness(2),throat_index),linspace(wallThickness(2),wallThickness(3),converge_index-throat_index),(wallThickness(3)*ones(1,(heightStepNumber - converge_index)))];

    
    %% Call Heat Transfer Solver
    [T_wgFinal, flowTemp, flowPressure, flowVelocity, finEfficiency, Qdot, finQdot, T_wl_Array, h_l_Array, h_g_Array] = SOLVER(chamberDiameterArray, heightStepNumber, newFluidProperties, channelWidths, channelHeights, wallThicknesses,inputValues);
    
    %% Structural Analysis
    chamberPressure = newFluidProperties(:,9)*100000; % Chamber Pressure (Pa)
    deltaP = flowPressure - chamberPressure';
    [vonMises,sigma_long, sigma_circ, sigma_rad, structuresOutput] = Structures(T_wl_Array, T_wgFinal, wallThicknesses, channelWidths, chamberPressure', chamberDiameterArray', deltaP, Qdot, inputValues);
    

    %% Plotter: Structures
    
    % Circ stress
        axialDist = newFluidProperties(:,1);
    
    
        figure;
        hold on; box on;
        
        % Plot with contrasting, muted colors (blue, orange, gray)
        p1 = plot(axialDist, sigma_circ,  '-', 'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]); % blue
        p2 = plot(axialDist, structuresOutput(2,:), '-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]); % orange
        p3 = plot(axialDist, structuresOutput(1,:),  '-', 'LineWidth', 1.8, 'Color', [0.50 0.50 0.50]); % gray
        
        % Axis labels
        xlabel('Location [m]', 'FontSize', 12);
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
        xlabel('Location [m]', 'FontSize', 12);
        ylabel('Channel Dimensions [mm]', 'FontSize', 12);
        
        % Grid and ticks
        grid on;
        grid minor;
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
                 'FontSize', 11, 'LineWidth', 1);
        xlim([0,axialDist(1)]);
        
        ylim([0, max(channelWidths)*1.1]);
        % Legend
        legend([h1 h2 h3], {'Channel Height', 'Channel Width', 'Wall Thickness'}, ...
               'Location', 'best', 'Box', 'off');
    
    % Von Mises vs Yield
        yieldStrength = (-(1.57*10^-5)*(T_wgFinal.^3) + (0.0139 * (T_wgFinal.^2)) - 4.21*(T_wgFinal) + 762)*1000000;
        UTS =  (-(1.31*10^-5)*(T_wgFinal.^3) + (0.0113 * (T_wgFinal.^2)) - 3.45*(T_wgFinal) + 719)*1000000;
            
        % Create figure
        figure;
        box on;
        
        % LEFT Y-AXIS: Strength & Stress [MPa]
        yyaxis left
        hold on;
        
        % Muted, contrasting colors: blue, orange, gray
        h1 = plot(axialDist, vonMises,'-', 'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]); % blue
        h2 = plot(axialDist, yieldStrength,'-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]); % orange
        h3 = plot(axialDist, UTS, '-', 'LineWidth', 1.8, 'Color', [0.50 0.50 0.50]); % gray
        
        ylabel('Strength & Stress [MPa]', 'FontSize', 12);
        
        
        allY_left = [UTS,vonMises];
        ylim([0, max(allY_left)*1.1]);
        
        
        %% RIGHT Y-AXIS: Contour (Radial Distance [m])
        yyaxis right
        
        h4 = plot(axialDist, chamberDiameterArray/2, 'k--', 'LineWidth', 1.8);
        ylabel('Radial Distance [m]', 'FontSize', 12);
        
        % --- Control RIGHT Y limits: 0 to next grid mark above max Contour ---
        
      
        ylim([0, max(chamberDiameterArray)*1.1]);
       
        
        %% Common X-axis
        xlabel('Axial Distance', 'FontSize', 12);
        
        % X limits: start and end at data bounds
        xlim([0,axialDist(1)]);
        
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


end