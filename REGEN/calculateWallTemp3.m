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
    [vonMises,sigma_long, sigma_circ, sigma_rad] = Structures(T_wl_Array, T_wgFinal, wallThicknesses, channelWidths, chamberPressure', chamberDiameterArray', deltaP, inputValues);


end