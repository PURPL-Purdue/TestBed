%% Main
function [] = CEAOut(p_c, OF, contour_name)
    
    % Get the folder where this script lives
    repoDir = fileparts(mfilename('fullpath'));

    % Add the subfolder dynamically
    addpath(fullfile(repoDir, 'MatlabCEA', '+CEA', 'bin'));
    
    recycle on;
    % Delete("CEAOut.xlsx");

    % Engine Variables
    % engineContour = readmatrix(contour_name);
    engineContour = readmatrix(contour_name);
    fuelTemperature = 293.15; %K
    oxTemperature = 293.15; %K
    subSonicARatios = engineContour(1:696,3)';%[linspace(10.496,1.1,300),linspace(1.101,1,200)];
    superSonicARatios = engineContour(697:end,3)';% [linspace(1.001,1.1,200),linspace(1.101,5.005,300)];
    compressionRatio = 11; %Unitless (for Finite Area Combustor option)
    mode = 'fz'; % Options are 'eq' and 'fz'
    FAC = 'false';

    %Engine Reactant Definitions
    reactants =   [                                         ...
                CEA.Reactant('RP-1',                        ...
                        'Type','Fuel',                     ...
                        'T',DimVar(fuelTemperature,'K'),   ...         
                        'Q',DimVar(1,'kg'))                ...
                CEA.Reactant('O2',                       ...
                        'Type','ox',                        ...
                        'T',DimVar(oxTemperature,'K'),      ...         
                        'Q',DimVar(1,'kg'))                 ...
                ];

    % CEA.Reactant('C3H8O,2-propanol',                           ...
    %                    'Type','Fuel',                     ...
    %                   'T',DimVar(fuelTemperature,'K'),             ...         
    %                  'Q',DimVar(1,'kg'),'E',-DimVar(200,'J/mol'),'rho',DimVar(0.786,'g/cm3'))

    a = 1;
    tempSub = [];
    tempSup = [];
    combined = [subSonicARatios,superSonicARatios];
    gigaMatrix = [];
    while a <= length(combined)/10
        if a <= length(subSonicARatios)/10
            tempSub = subSonicARatios(((a-1)*10)+1:a*10);
            tempSup = [];
        else
            tempSup = combined((((a-1)*10)+1:a*10));
            tempSub = [];
        end
        if a >= length(subSonicARatios)/10 && mode == "eq"
            mode = 'fz';
        end
    %% CEA call
        if mode == "eq"
            data =  CEA.Run(reactants,                                ...
                'ProblemType','Rocket',                         ...
                'Flow',mode,'CombLength','fac'  ,                               ...
                'Pc',DimVar(p_c,'psi'),                           ...
                'OF',OF,'CR',compressionRatio, 'subar',tempSub,'supar',tempSup,                                       ...
                'Outputs',{'mach','gamma','t','visc','Cp','p','Ae/At','Prandtl'});
            FAC = 'true';
        elseif mode == "fz"
            data =  CEA.Run(reactants,                                ...
                'ProblemType','Rocket',                         ...
                'Flow',mode ,                               ...
                'Pc',DimVar(p_c,'psi'),                           ...
                'OF',OF,'CR',compressionRatio, 'subar',tempSub,'supar',tempSup,                                       ...
                'Outputs',{'mach','gamma','t','visc','Cp','p','Ae/At','Prandtl'});
        else
            fprintf("Mode error")
            break;
        end
        
        %% Data extraction
        rawData = readtable("Detn.out", 'FileType','text', 'VariableNamingRule','preserve');
        
        %Variables to extract
        targets = ["Ae","PRANDTL","MACH","GAMMAs","T,","VISC,MILLIPOISE","Cp,","P,","CSTAR,","RHO,","SON","Isp,"];
        
        % Data reader
        out = [];
        outCol = 1;
        
        targetNum = 1;
        while targetNum <= length(targets)

            line = 1;
            while line < height(rawData(:,1))
                currentLine = split(string(rawData{line,1})," ");
                if currentLine(1) == targets(targetNum)
                    rowNum = 1;
                    %Number extraction
                    if string(rawData{line,2}) == ""
                        rawValueLine = split(string(rawData{line,1})," ");
                        rawValueLine = rawValueLine(3:end);
                        valueLine = [];
                        i = 1;
                        while i <= length(rawValueLine)
                            if rawValueLine(i) ~= ""
                                valueLine(end+1) = rawValueLine(i);
                            end
                            i = i + 1;
                        end
                    else
                        rawValueLine = split(string(rawData{line,2})," ");
                        rawValueLine = rawValueLine(3:end);
                        valueLine = strings(1,1);
                        i = 1;
                        while i <= length(rawValueLine)
                            if rawValueLine(i) ~= "" && rawValueLine(i) ~= "0"
                                valueLine(end+1) = rawValueLine(i);
                            end
                            i = i + 1;
                        end
                    end
                    %Removes injector values
                    while length(valueLine) > 10
                        valueLine = valueLine(2:end);
                    end               
                    %Searches for exponential notation and converts strings to
                    %usable numbers
                    q = 1;
                    while q <= length(valueLine)
                        if isstring(valueLine(q)) && contains(valueLine(q),"-") == true
                            notation = 10 .^ -double((extractBetween(valueLine(q),strfind(valueLine(q),"-")+1,length(valueLine)-2)));
                            newDouble = notation * double(extractBetween(valueLine(q),1,strfind(valueLine(q),"-")-1));
                            valueLine(q) = newDouble;
                        end
                        q = q + 1;
                    end
                    %Add to out
                    while rowNum <= length(valueLine)
                        out(rowNum,outCol) = valueLine(rowNum);
                        rowNum = rowNum + 1;
                    end
                    outCol = outCol + 1;
                    break;
                end
                    line = line + 1;
            end
            targetNum = targetNum + 1;
        end
        %Save output to temperary array
        gigaMatrix = vertcat(gigaMatrix, out);
        clc;
        fprintf(a + "/100 CEA iterations finished")
        a = a + 1;
    end

    %% Add Axial positions to area ratios
    contour = [engineContour(:,4:5), engineContour(:,3)];
    %gigaMatrix = [gigaMatrix(1:695,:);gigaMatrix(708:end,:)];
    %contour = readmatrix("Engine Contour Cleaned and Sorted (Metric).csv");

    %gigaMatrix = [gigaMatrix(1:695,:);gigaMatrix(708:end,:)];

    subSonicContour = contour(1:696,:); % contour(1:783,:);
    supSonicContour = contour(697:end,:); % contour(784:end,:);

    axPos = [];
    change = 0;
    currentOutNum = 1;
    while currentOutNum <= length(gigaMatrix(:,1))
        currentARatio = gigaMatrix(currentOutNum,1);
        currentAxNum = 2;
        if change == 0
            while currentAxNum <= length(subSonicContour(:,1))
                if subSonicContour(currentAxNum,3) <= currentARatio
                    %Linear interpoliation
                    delA = subSonicContour(currentAxNum,3)-subSonicContour(currentAxNum-1,3);
                    delAx = subSonicContour(currentAxNum,1)-subSonicContour(currentAxNum-1,1);
                    if delA == 0
                        delA = 1;
                        delAx = 0;
                    end
                    axPos(end+1) = (delAx/delA) * (currentARatio-subSonicContour(currentAxNum-1,3)) + subSonicContour(currentAxNum,1);
                    break;
                end
                currentAxNum = currentAxNum + 1;
            end
        else
            while currentAxNum <= length(supSonicContour(:,1))
                if supSonicContour(currentAxNum,3) >= currentARatio
                    %Linear interpolation
                    delA = supSonicContour(currentAxNum,3)-supSonicContour(currentAxNum-1,3);
                    delAx = supSonicContour(currentAxNum,1)-supSonicContour(currentAxNum-1,1);
                    if delA == 0
                        delA = 1;
                        delAx = 0;
                    end
                    axPos(end+1) = (delAx/delA) * (currentARatio-supSonicContour(currentAxNum-1,3)) + supSonicContour(currentAxNum,1);
                    break;
                end
                currentAxNum = currentAxNum + 1;
            end
        end
        if currentARatio == 1
            change = 1;
        end
        currentOutNum = currentOutNum + 1;
    end

    gigaMatrix = [contour(:,1), gigaMatrix];

    %Write output to table
    temp = [0,"Ae/At","PRANDTL","MACH","GAMMAs","T,","VISC,MILLIPOISE","Cp,","P,","CSTAR,","RHO,","SON","Isp,"];
    options = ["",sprintf("Mode: %s",mode),sprintf("Fac: %s",FAC),sprintf("Compression Ratio: %f",compressionRatio),sprintf("OF: %f",OF),sprintf("P_c: %f",p_c),sprintf("Point Density: 1000"),"","","","","",""];
    writematrix([options;temp;gigaMatrix],"CEA_Maelstrom.xlsx")

    fprintf("CEA Values Obtained Successfully.")

end