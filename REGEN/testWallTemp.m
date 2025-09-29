% PSP dimensions
W1 = 0.1772*0.0254;
H1 = 0.1181*0.0254;
T1 = 0.0394*0.0254;
W2 = 0.1772*0.0254;
H2 = 0.1181*0.0254;
T2 = 0.0394*0.0254;
W3 = 0.0591*0.0254;
H3 = 0.0787*0.0254;
T3 = 0.0295*0.0254;
W4 = 0.0894*0.0254;
H4 = 0.0894*0.0254;
T4 = 0.0394*0.0254;

combo1 = [W1, H1, T1];
combo2 = [W2, H2, T2];
combo3 = [W3, H3, T3];
combo4 = [W4, H4, T4];

heightStepNumber = 67;
chamberLength = 8.97*0.0254;

heightStepArray = chamberLength/heightStepNumber:chamberLength/heightStepNumber:chamberLength;
channelMatrix = [];

throatInd = 1+ length(heightStepArray(heightStepArray < 0.178816));
startConvInd = 1+length(heightStepArray(heightStepArray < 0.132334));


channelMatrix(startConvInd:throatInd,1) = linspace(W3,W2,throatInd-startConvInd+1);
channelMatrix(startConvInd:throatInd,2) = linspace(H3,H2,throatInd-startConvInd+1);
channelMatrix(startConvInd:throatInd,3) = linspace(T3,T2,throatInd-startConvInd+1);

channelMatrix(throatInd:heightStepNumber,1) = linspace(W4,W3,heightStepNumber-throatInd+1);
channelMatrix(throatInd:heightStepNumber,2) = linspace(H4,H3,heightStepNumber-throatInd+1);
channelMatrix(throatInd:heightStepNumber,3) = linspace(T4,T3,heightStepNumber-throatInd+1);
i = 1;

for a = heightStepArray
    
    if(heightStepArray(i)<=startConvInd)
        channelMatrix(i,4) = 0.09525; % set diameter in m
        channelMatrix(i,1:3) = combo1;
    elseif(heightStepArray(i)<=throatInd)
        channelMatrix(i,4) = (((heightStepArray(i)-0.132334)*-1.245) +0.09525);
    else
        channelMatrix(i,4) = (((heightStepArray(i)-0.178816)*0.5129) +0.037338);
    end

    i=i+1;
end



numchannels = 62;


%% Inlet Condition Values
T_start= 292; % K
P_start = 2517000; % Pa
v_start = 2.71; %m/s calculated from mass flow through channels at start
rho_start = 999.63; %kg/m^3
m_flow_total = 2.26796; %kg/s
mass_flow = m_flow_total/numchannels;

%% Engine values
chamberPressure = 1379000; % Chamber Pressure (Pa)
k_w = 162; % thermal conductivity of the wall (W/m*K)
gravity = 9.81; %m/s^2

%% Array initialization
% initialize/re-define matrices & arrays
flowTemp = [];
flowVel = [];
flowPressure = [];
currentHeightStep = heightStepArray(end)/heightStepNumber; % axial length of current height step
targetTemp = 573 ; % K temp of 50%yield strength decrease in al6061 ram2 (https://www.elementum3d.com/wp-content/uploads/2024/02/A6061-RAM2-2-pg-Web-Event-Data-Sheets-2024-01-18.pdf). Setting max hotwall-gas side temp
heightStepNumber = 0; % initialize height step number

fluidProperties = readmatrix("CEAOutFz_PSP.xlsx"); %pull all nasaCEA values into fluidProperties
    fluidProperties(1,:) = [];
    y = 1;
    r = 1;
    axialDist = (fluidProperties(:,1));
    newFluidProperties = zeros(length(heightStepArray),10);
    newFluidProperties(:,1) = heightStepArray;
    while y <= length(heightStepArray) % translating CEA outputs to height step number length output by averaging values over height step number
        
        
        a = r; % MAY NEED TO CHANGE BASED ON WHAT GETS READ FROM EXCEL FILE (add 2 or something to accoutn for text)

        sumAEAT = 0;
        sumPrandtl = 0;
        sumMach = 0;
        sumGamma = 0;
        sumT = 0;
        sumVisc = 0;
        sumCp = 0;
        sumP = 0;
        sumCstar = 0;                

        while a <= length(axialDist)
            if axialDist(a) < heightStepArray(y)
                sumAEAT = sumAEAT+fluidProperties(a,2);
                sumPrandtl = sumPrandtl+fluidProperties(a,3);
                sumMach = sumMach+fluidProperties(a,4);
                sumGamma = sumGamma+fluidProperties(a,5);
                sumT = sumT+fluidProperties(a,6);
                sumVisc = sumVisc+fluidProperties(a,7);
                sumCp = sumCp+fluidProperties(a,8);
                sumP = sumP+fluidProperties(a,9);
                sumCstar = sumCstar+fluidProperties(a,10);


                
                a=a+1;
            else
                divFactor = a-r;
                newFluidProperties(y,2) = sumAEAT/divFactor;
                newFluidProperties(y,3) = sumPrandtl/divFactor;
                newFluidProperties(y,4) = sumMach/divFactor;
                newFluidProperties(y,5) = sumGamma/divFactor;
                newFluidProperties(y,6) = sumT/divFactor;
                newFluidProperties(y,7) = sumVisc/divFactor;
                newFluidProperties(y,8) = sumCp/divFactor;
                newFluidProperties(y,9) = sumP/divFactor;
                newFluidProperties(y,10) = sumCstar/divFactor;
                r = a;
                
                break; 
            
            end
        end
        y=y+1;
    end

wallThicknesses = 0.001:0.001:0.02; %Prospective wall thicknesses range in meters
gasTemps = newFluidProperties(:,6); %Hot gas temperature in K
f = H_g_From_Temperature(targetTemp,newFluidProperties); % for debugging, delete when done  
flux = f .* (gasTemps - targetTemp); %Calcs heat flux in W/m^2

        
    %Outputs
    QDot = flux;
    

%% Main Height Step Loop
for i = heightStepArray
    heightStepNumber = heightStepNumber + 1; % iterate through each height step
    
    width = channelMatrix(heightStepNumber,1);
    height = channelMatrix(heightStepNumber,2);
    wall_thickness = channelMatrix(heightStepNumber,3);
    chamberDiameter = channelMatrix(heightStepNumber,4);
    Tc = -(wall_thickness * flux / k_w - (targetTemp)); %Calcs coolant wall temperature required to induce heat transfer
    Tsubc = Tc;
    if (heightStepNumber==1)

        hotWall_dP = P_start - chamberPressure; %calculate dP for structures (Pa)
        pressure = P_start;
        temp = T_start;
        velocity = v_start;
    else

        hotWall_dP = flowPressure(heightStepNumber-1) - chamberPressure; %calculate dP for structures (Pa)
        pressure = flowPressure(heightStepNumber-1);
        temp = flowTemp(heightStepNumber-1);
        velocity = flowVel(heightStepNumber-1);
    end
    %% Call Tucker's Function HERE, update variables (coolant side hotwall temp, Heat flux, Wall thickness) (needs updated wall and dP)
    bendMaxs = 3 * hotWall_dP* width ^ 2 / (4 * wall_thickness^ 2); %Calcs expected max bending stress in Pa
    TMaxs = 3 * hotWall_dP* width / (0.577 * 4 * wall_thickness); %Calcs expected tensile stress in Pa


    if(temp == 999999999999999999999999999999999999999999) % if channel dimension combo is already unsuccessful, do not let computation with it continue
        flowTemp(heightStepNumber) = 999999999999999999999999999999999999999999;
        break
    end

    %% Reynold's Number
    % Hydraulic Diameter (m)
    hyd_diam = (2*width*height)/(height+width);
    

    % Density (kg/m^3)
    density = py.CoolProp.CoolProp.PropsSI('D','P',pressure,'T',temp,'Water');

    % Dynamic viscosity [Pa·s]
    dyn_visc = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',pressure,'T',temp,'Water');%dynamic viscosity is normal viscosity, for kinematic viscosity, divide by density


    % Calculate Reynolds number
    Re = (density * velocity * hyd_diam) / dyn_visc;
    


    %% Calculate Prandtl Number
    % Thermal conductivity [W/(m·K)]
    kf = py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY','P',pressure,'T',temp,'Water');
    %^Use empirical data/curves found from papers in regen channel

    % Specific heat capacity [J/(kg·K)]
    cp = py.CoolProp.CoolProp.PropsSI('C','P',pressure,'T',temp,'Water');
    
    %^Use empirical data/curves found from papers in regen channel

    % Prandtl number [-]
    Pr = cp * dyn_visc / kf;



    %% Sieder Tate Nusselt's Number
    % Get viscosity at wall temperature for Sieder-Tate correction
    mu_wall = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',pressure,'T',Tsubc(heightStepNumber),'Water');
    %(1.74 + -0.0493*Tsubc(heightStepNumber) + 6.98*(10^-4)*(Tsubc(heightStepNumber)^2) -3.78*(10^-6)*(Tsubc(heightStepNumber)^3))/1000;

    % Calculate Nusselt number using Sieder-Tate correlation
    Nu = 0.027 * Re^(4/5) * Pr^(1/3) * (dyn_visc/mu_wall)^0.14;
   
    % Calculate convective heat transfer coefficient [W/(m^2·K)]
    h_l = Nu * kf / hyd_diam;
    


    %% Calculate Required flow temperature
    % Area of Fin (m^2)
    A_fin = currentHeightStep*2*height;

    % Area of Wall on Coolant side (m^2)
    A_wallL = currentHeightStep*width;

    % Calculate fin width at start of channels (m)
    fin_width = ((pi*(chamberDiameter+2*(wall_thickness))) - (numchannels*width))/numchannels;
    
    % Area of Wall on Hotwall side (m^2)
    A_wallG = currentHeightStep *(width+fin_width);
    
    % Finning Parameter- measure of convection from fins
    fin_param = sqrt((2*h_l*(fin_width+currentHeightStep))/(k_w*currentHeightStep*fin_width));

    % Fin Efficiency- fin convection efficiency
    fin_efficiency = tanh(fin_param*height)/(fin_param*height);

    % Calculate angle of channel slice
    angle_channel = (width/(pi*(chamberDiameter+2*(wall_thickness))))*360;

    % Calculate mass of coolant in height step, m_coolant (kg)
    m_coolant = density*currentHeightStep*((pi*((chamberDiameter+wall_thickness*2)/2)^2)-(pi*((chamberDiameter/2)^2)))*(angle_channel/360);

    % Calculate required specific heat transfer rate, qdotL_total (J/kg*s)
    qdotL_total = QDot(heightStepNumber)*((pi*(chamberDiameter)*(angle_channel/360))*currentHeightStep)/m_coolant;

    % Calculate required coolant temp, T_L_req (K)
    T_L_req =  -((qdotL_total*A_wallG)/(h_l*((fin_efficiency*A_fin)+A_wallL)))+ Tsubc(heightStepNumber);



    %% Calculate Coolant Temp increase, delta_T (K)
    delta_T = qdotL_total*currentHeightStep*(width+(2*height*fin_efficiency))/(mass_flow*cp);

    if (heightStepNumber == 1)
        flowTemp(heightStepNumber) = T_start + delta_T; % initial flow temp plus deltaT

    else
        flowTemp(heightStepNumber) = flowTemp(heightStepNumber-1)+delta_T; % Previous step flow temp plus deltaT

    end

    if (temp > T_L_req)
        flowTemp(heightStepNumber) = 999999999999999999999999999999999999999999; % if coolant temp is too high, nullify

    end


    %% Calculate Coolant Pressure Drop
    if(Re < 2100) % Calculate Friction Factor, cf
        cf = 16/Re;
    elseif(5000 < Re && 200000 > Re)
        cf = 0.046/(Re^0.2);
    else
        cf = 0.0014 + (0.125/(Re^0.32));

    end

    delta_P = 2 * cf * currentHeightStep * density * (velocity^2)/ hyd_diam; %Frictional static pressure drop across the channel
    
    if(heightStepNumber ==1)
        flowPressure(heightStepNumber) = P_start;
    else
        flowPressure(heightStepNumber) = flowPressure(heightStepNumber-1) - delta_P; %Flow pressure
    end 
    %% Calculate Coolant Velocity Increase via Bernoulli's
    flowVel(heightStepNumber) = sqrt((delta_P/density)+(gravity*(currentHeightStep))+velocity^2); % Flow velocity
    display(T_L_req)
end

%display(flowTemp)
%display(flowVel)
%display(flowPressure)