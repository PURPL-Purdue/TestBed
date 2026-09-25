function [] = generateContour(hotwallGeometry, filename, resolution)
    
    in_to_m = 0.0254;

    % Extract parameters
    chamberDiameter = hotwallGeometry(1);
    throatDiameter  = hotwallGeometry(2);
    exitDiameter    = hotwallGeometry(3);
    convergingAngle = hotwallGeometry(4);
    divergingAngle  = hotwallGeometry(5);
    totalLength     = hotwallGeometry(6);
    convergingFillet= hotwallGeometry(7);
    throatFillet    = hotwallGeometry(8);

    % Convert to meters
    cR = chamberDiameter/2 * in_to_m;
    tR = throatDiameter/2   * in_to_m;
    eR = exitDiameter/2     * in_to_m;
    tL = totalLength        * in_to_m;
    cF = convergingFillet   * in_to_m;
    tF = throatFillet       * in_to_m;

    % Convert angles
    cA = deg2rad(convergingAngle);
    dA = deg2rad(divergingAngle);

    % Generate points
    x_points = linspace(0, tL, resolution);
    r_points = zeros(1, resolution);

    for i = 1:resolution
        r_points(i) = contour(x_points(i), cR, tR, eR, tL, cF, tF, cA, dA);
    end

    %──────────────────────────────────────────────────────────────
    % FIND THROAT RADIUS DIRECTLY FROM CONTOUR DATA
    %──────────────────────────────────────────────────────────────
    r_throat = min(r_points);      % NO THROAT DIAMETER USED
    area_ratio = (r_points ./ r_throat).^2;

    %──────────────────────────────────────────────────────────────
    % EXPORT: x, r, area_ratio (NO HEADERS, RAWDOGGED)
    %──────────────────────────────────────────────────────────────
    outData = [x_points(:), r_points(:), area_ratio(:)];

    repoDir = fileparts(mfilename('fullpath'));   % folder the .m file lives in
    outName = fullfile(repoDir, "Contour_" + filename + ".xlsx");

    % Delete old file if it exists
    if isfile(outName)
        delete(outName);
    end

    writematrix(outData, outName, 'FileType', 'spreadsheet');


    %──────────────────────────────────────────────────────────────
    % PLOT
    %──────────────────────────────────────────────────────────────
    % figure;
    % plot(x_points, r_points, 'b-', 'LineWidth', 2);
    % xlabel('Axial Distance (m)');
    % ylabel('Radial Distance (m)');
    % title('Nozzle Contour');
    % grid on;
    % axis equal;
end



%% ─────────────────────────────────────────────────────────────
%  CONTOUR FUNCTION
% ─────────────────────────────────────────────────────────────
function r = contour(x, cR, tR, eR, tL, cF, tF, cA, dA)

    h  = (tF + tR - tF*cos(dA) - eR + tL*tan(dA))/tan(dA) - tF*sin(dA);
    k  = tR + tF;

    xc = h - tF*sin(cA) ...
         + (cR - (k - tF*cos(cA)))/(-tan(cA)) ...
         + cF/tan((pi + cA)/2);

    rc = cR - cF;

    A = k - tF*cos(cA);
    m = tan(cA);
    B = h - tF*sin(cA);

    E = A - rc + m*B;

    a = m^2 + 1;
    b = -(2*xc + 2*E*m);
    c = (E^2 - cF^2 + xc^2);

    x_int1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
    x_int2 = (-b - sqrt(b^2 - 4*a*c)) / (2*a);
    x_filletEnd = min(x_int1, x_int2);

    x_flatEnd      = xc;
    x_convergeStart= B;
    x_throat       = (tF + tR - tF*cos(dA) - eR + tL*tan(dA))/tan(dA);


    if x < x_flatEnd
        r = cR;
    elseif x < x_filletEnd
        r = rc + sqrt(cF^2 - (x - xc).^2);

    elseif x < x_convergeStart
        r = A - m*(x - B);

    elseif x < x_throat
        r = k - sqrt(tF^2 - (x - h).^2);

    elseif x <= tL
        r = tan(dA)*x + (eR - tL*tan(dA));
    end
end