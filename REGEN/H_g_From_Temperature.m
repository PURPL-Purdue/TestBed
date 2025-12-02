function h_out = H_g_From_Temperature(temperature, fluid_information, inputFlowValues)
    

    info = fluid_information';
    T_stag = info(6,:) .* (1+((info(5,:)-1)/2).*(info(4,:).^2));
    
    %% CRIIIIINGE IMPERIAL :(    
    %Input manipulation
    wall_temp = temperature * 1.8; %Conversion to Rankine
    
    R = inputFlowValues(12); %inches
    G = 32.17405; %ft/s^2
    throatDiameter = inputFlowValues(11); %inches


    area_Ratio = info(2,:);

    %Unit conversion
    info(6,:) = (T_stag.* 1.8); %Conversion to Rankine
    info(7,:) = info(7,:) .* 0.0000671968994813 ./ 12; %Conversion to lb/in-sec 

    info(8,:) = info(8,:) .* 0.2388458966; %Conversion to BTU/lb-F
    info(9,:) = info(9,:) .* 14.503773800722; %Conversion to psia
    info(10,:) = info(10,:) .* 3.2808399; %Conversion to ft/s


    %Property Variation calc
    prop = 1 ./ ((0.5 .* ((wall_temp)./(info(6,:))) .* (1 + ((info(5,:)-1) ./ 2) .* info(4,:))+0.5).^ 0.68 .* (1 + (info(5,:) - 1)./(2) .* info(4,:)).^0.12);

    %H calc
    hsubg = ((0.026./((throatDiameter.^0.2))) .* ((info(7,:).^0.2 .* info(8,:))./(info(3,:).^0.6)) .* ((info(9,:) .* G) ./ (info(10,:))).^0.8 .* (throatDiameter./(R)).^0.1) .* (info(2,:).^(-1)).^0.9 .* prop;

    %Output manipulation
    hsubg = hsubg .* 2941643.128651; %Btu/s * in * F to Watt/M*K = 1055.06 *39.37 * 9/5

                                     % 


    h_out = hsubg';

    %% Based Metric :)
    % %Input manipulation
    % wall_temp = temperature; %Conversion to Rankine
    % info = fluid_information';
    % R = 0.61625/39.37; %m
    % G = 9.81; %m/s^2
    % throatDiameter = 0.92596/39.37; %m
    % 
    % 
    % area_Ratio = info(2,:);
    % 
    % %Unit conversion
    % 
    % info(8,:) = info(8,:) .* 1000; %Conversion to J/kg*K
    % 
    % %Property Variation calc
    % prop = 1 ./ ((0.5 .* ((wall_temp)./(info(6,:))) .* (1 + ((info(5,:)-1) ./ 2) .* info(4,:))+0.5).^ 0.68 .* (1 + (info(5,:) - 1)./(2) .* info(4,:)).^0.12);
    % 
    % %H calc
    % hsubg = ((0.026./((throatDiameter.^0.2))) .* ((info(7,:).^0.2 .* info(8,:))./(info(3,:).^0.6)) .* ((info(9,:)) ./ (info(10,:))).^0.8 .* (throatDiameter./(R)).^0.1) .* (info(2,:).^(-1)).^0.9 .* prop;
    % 
    % % hsubg = (0.026./((throatDiameter.^0.2))) .* ((info(7,:).^0.2 .* info(8,:))./(info(3,:).^0.6)) .* (((info(9,:)) ./ (info(10,:))).^0.8) .* ((throatDiameter./R).^0.1) .* ((info(2,:).^(-1)).^0.9) .* prop;
    % 
    % %(((info(7,:).^0.2) .* info(8,:))./((info(3,:).^0.6)))
    % 
    % 
    % 
    % h_out = hsubg';
end
