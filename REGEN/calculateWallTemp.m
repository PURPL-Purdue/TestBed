function [flowTemp, flowVel] = calculateWallTemp(heightStepNumber, heatFlux, wallTempArray, hotwallSurfaceAreaArray, flowTempMatrix, heightArray, widthArray)
    
    flowTemp = zeros([length(widthArray), length(heightArray)]); %Preallocate array sizes
    flowVel = flowTemp; %Preallocate array sizes
    MASS_FLOW; % Precalcuated mass flow based on # of channels in Malestrum
    wInd = 0;
    hInd = 0;
    
    for width = widthArray
        wInd = wInd + 1;
        for height = heightArray
            hInd = hInd + 1;
            flowVel(wInd, hIn) = 1; % Flow velocity
            flowTemp(wInd, hIn) = 1; %Flow temp
            coolentMass; %%Calculated from volume and denisty
            fluidConactArea;
            convtCoeff;
        end
    end

end