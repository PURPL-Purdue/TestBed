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
        xlabel('Location [in]', 'FontSize', 12);
        ylabel('Circumferential stress [MPa]', 'FontSize', 12);
        
        % Grid and ticks
        grid on;
        grid minor;
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
                 'FontSize', 11, 'LineWidth', 1);
        xlim([axialDist(1) axialDist(end)]);
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
        xlim([axialDist(1) axialDist(end)]);
        % Legend
        legend([h1 h2 h3], {'Channel Height', 'Channel Width', 'Wall Thickness'}, ...
               'Location', 'best', 'Box', 'off');
    



end