%% Main
%function [out] = CEAOut()
clc;
clear;

recycle on;
%delete("CEAOut.xlsx");

%Engine Variables
engineContour = readmatrix("PSP_Engine_Contour.xlsx");
%engineContour = readmatrix("Engine Contour Cleaned and Sorted (Metric).csv");
p_c = 250; %PSI
OF = 1.2; %Unitless
fuelTemperature = 293.15; %K
oxTemperature = 90; %K
subSonicARatios = engineContour(1:790,3)';%[linspace(10.496,1.1,300),linspace(1.101,1,200)];
superSonicARatios = engineContour(779:end,3)';% [linspace(1.001,1.1,200),linspace(1.101,5.005,300)];
compressionRatio = 11; %Unitless (for Finite Area Combustor option)
throatArea = 0.00097314; %m^2
mode = 'fz'; %Options are 'eq' and 'fz'
FAC = 'false';


%Engine Reactant Definitions
reactants =   [                                           ...
            
%CEA.Reactant('C3H8O,2-propanol',                           ...
 %                   'Type','Fuel',                     ...
  %                   'T',DimVar(fuelTemperature,'K'),             ...         
   %                  'Q',DimVar(1,'kg'),'E',-DimVar(200,'J/mol'),'rho',DimVar(0.786,'g/cm3'))


CEA.Reactant('C2H5OH(L)',                           ...
                     'Type','Fuel',                     ...
                     'T',DimVar(fuelTemperature,'K'),             ...         
                     'Q',DimVar(1,'kg'),'E',-DimVar(1367,'kJ/mol'))                ...
            CEA.Reactant('O2(L)',                            ...
                    'Type','ox',                        ...
                    'T',DimVar(oxTemperature,'K'),              ...         
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
    rawData = readtable("Detn.out",'Filetype','text');
    
    
    %Variables to extract
%    targets = ["Ae/At",,"PRANDTL","MACH","GAMMAs","T","VISC","Cp","P","CSTAR","RHO","SON","Isp"];
    targets = ["Ae/At","Ae","PRANDTL","MACH","GAMMAs","T,","VISC,MILLIPOISE","Cp,","P,","CSTAR,","RHO,","SON","Isp,"];
    
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
    a = a + 1;
end

%% Add Axial positions to area ratios
contour = [engineContour(:,4:5), engineContour(:,3)];
%gigaMatrix = [gigaMatrix(1:695,:);gigaMatrix(708:end,:)];
%contour = readmatrix("Engine Contour Cleaned and Sorted (Metric).csv");

 %gigaMatrix = [gigaMatrix(1:695,:);gigaMatrix(708:end,:)];

subSonicContour = contour(1:784,:)%contour(1:783,:);
supSonicContour = contour(785:end,:)%contour(784:end,:);

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
gigaMatrix = [[contour(1:end-2,1)']',gigaMatrix];

%Write output to table
temp = [0,targets];
options = ["",sprintf("Mode: %s",mode),sprintf("Fac: %s",FAC),sprintf("Compression Ratio: %f",compressionRatio),sprintf("OF: %f",OF),sprintf("P_c: %f",p_c),sprintf("Point Density: 1000"),"","","","","",""];
writematrix([options;temp;gigaMatrix],"CEAOutFzPsp_IAC_11-23-25.xlsx")
%end