function h_out = H_g_From_Temperature(temperature, fluid_information)

    %Input manipulation
    wall_temp = temperature * 1.8 + 491.67; %Conversion to Rankine
    info = fluid_information';
    
    area_Ratio = info(2);
    
    %Unit conversion
    info(6) = (info(6).* 1.8); %Conversion to Rankine
    info(8) = info(8) .* 0.0000671968994813 ./ 12; %Conversion to lb/in-sec
    info(9) = info(9) .* 0.2388458966; %Conversion to BTU/lb-F
    info(12) = info(12) .* 3.2808399; %Conversion to ft/s
    info(10) = info(10) .* 14.503773800722; %Conversion to psia
    
    %Property Variation calc
    prop = 1 ./ ((0.5 .* ((wall_temp)./(info(6))) .* (1 + ((info(5)-1) ./ 2) .* info(4))+0.5).^ 0.68 .* (1 + (info(5) - 1)./(2) .* info(4)).^0.12);
        
    %H calc
    hsubg = ((0.026./((info(7).^0.2))) .* ((info(8).^0.2 .* info(9))./(info(3).^0.6)) .* ((info(10) .* info(11)) ./ (info(12))).^0.8 .* ((info(7))./(info(13))).^0.1) .* (info(2).^(-1)).^0.9 .* prop;
    
    %Output manipulation
    hsubg = hsubg .* 2941643.128651;
    
    
    h_out = hsubg;
end
