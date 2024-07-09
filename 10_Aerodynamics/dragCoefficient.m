% Calculates drag coeffcient and angle of attack from given C_L in AVL
% Inputs: airspeed, lift coefficient (both can be an array, but have to have the same size!)
% Outputs: Drag Coefficient and Angle of Attack
% TO DO:
function [C_D, AoA] = dragCoefficient(aircraftdata,inputdata,airspeed,C_L)

if any(strcmp(fieldnames(inputdata),'offsetdata'))
    offsetdata = inputdata.offsetdata;
    if any(strcmp(fieldnames(offsetdata),'drag'))
        offsetFactor = offsetdata.drag/100 + 1;
    else
        offsetFactor = 1;
    end
else
    offsetFactor = 1;
end

%% Extract data from aircraftdata
b_wing = aircraftdata.Geometry.wing.span;
taper_wing = aircraftdata.Geometry.wing.taper_ratio;
S_ref = aircraftdata.Geometry.wing.projected_surface;
airfoilName = aircraftdata.Geometry.wing.airfoil;
C_D_other = aircraftdata.Aerodynamics.drag_coefficient_other;

%% Extract data from inputdata
viscosity = inputdata.viscosity;
density = inputdata.density;
Airfoil_Data_Struct = inputdata.airfoil_data.(airfoilName).polar;

C_D = zeros(1,length(C_L));
AoA = zeros(1,length(C_L));

if length(airspeed) ~= length(C_L)
    if length(airspeed) == 1
        airspeed = ones(1,length(C_L)) * airspeed;
    elseif length(C_L) == 1
        C_L = ones(1,length(airspeed)) * C_L;
    else
        error('Airspeed array and C_L array have to have the same length')
    end
end
%% Calculate drag for all C_L's and airspeeds
for i = 1:length(C_L)
    if airspeed(i) <= 0
        C_D(i) = 0;
        continue
    end
    % Calculate induced drag coefficient, if C_L is not 0
    if C_L(i) ~= 0
        [y,C_l,~,~,C_D_i,AoA(i),stall] = AVL_Call(aircraftdata,inputdata,airspeed(i),C_L(i));
    else
        y = zeros(1,2);
        C_l = zeros(1,2);
        C_D_i = 0;
        stall = false;
    end
    % Calculate chord length for various y-coordinates
    %y = linspace(0,b_wing/2,10);
    chord = 2 * S_ref / (1 + taper_wing) / b_wing * (1 - (1 - taper_wing)/b_wing*2*y);
    % Calculate Zero-Lift drag coefficient for all sections of the wing
    C_D_0 = zeros(1,length(y)-1);
    for j = 1:length(y)-1
        average_chord = (chord(j)+chord(j+1))/2;
        average_Cl = (C_l(j)+C_l(j+1))/2;
        area = 2 * average_chord * (y(j+1)-y(j));
        Re = density * airspeed(i) * average_chord / viscosity;
        if Re < 5000
            C_D_0_section = 0;
        else
            [C_D_0_section, ~] = Cd_Cm_of_Cl(Airfoil_Data_Struct, 0, Re);
        end
        C_D_0(j) = C_D_0_section * area / S_ref;
    end
    C_D_0 = sum(C_D_0);
    if stall
        error('Error105:dragStall',['Stall occured when calculating drag at airspeed of ',num2str(airspeed(i)),' m/s and a C_L of ',num2str(C_L(i))]);
    end
    if isnan(C_D_other)
        C_D_other = 0;
    end
    C_D(i) = (C_D_i + C_D_0 + C_D_other) * offsetFactor;
end
end
