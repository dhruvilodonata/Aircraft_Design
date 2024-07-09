% Calculates the lift distribution for max load case and other cases
% Calculates maximum achievable lift coeffcient at max Speed
% TO DO:
% 
function [aircraftdata,error_cell] = liftDistribution(aircraftdata,inputdata,error_cell)
%% Import of Input Data
density = inputdata.density;
g = inputdata.gravity;

%% Extract data drom aircraft data struct
speed_stall = aircraftdata.Performance.stall_speed;
speed_cruise = aircraftdata.Performance.cruise_speed;
speed_pattern = inputdata.speed_pattern;
speed_maxLoad = aircraftdata.Performance.maximum_speed;

S_ref = aircraftdata.Geometry.wing.projected_surface;
mass = aircraftdata.Mass.total_mass;

%% Calculation of max C_L
C_L_max = max3DLiftCoefficient(aircraftdata, inputdata, speed_maxLoad);

%% Generation of lift distribution
cases = fields(aircraftdata.Aerodynamics.wing_loading);
if isnan(speed_stall) || isnan(speed_cruise)
    speeds = [speed_pattern, speed_maxLoad];
    cases = cases(3:4,1);
else
    speeds = [speed_stall, speed_cruise, speed_pattern, speed_maxLoad];
end

C_Ls = 2 * mass * g ./ (density * speeds.^2 * S_ref);
C_Ls(end) = C_L_max;


for i = 1:length(speeds)
    caseName = cases{i};
    V = speeds(i);
    C_L = C_Ls(i);
    [y,C_l,C_m,l_C_l,~,~,stall] = AVL_Call(aircraftdata,inputdata,V, C_L);
    if stall
        error('Error105:stall',['Stall occured when calculating lift distribution at airspeed of ',num2str(V),' m/s and a C_L of ',num2str(C_L)]);
    end
    for j = 1:length(y)
        aircraftdata.Aerodynamics.wing_loading.(caseName).load_table(j,:) = {y(j), C_l(j), C_m(j), l_C_l(j)};
    end
    aircraftdata.Aerodynamics.wing_loading.(caseName).speed = speeds(i);
    aircraftdata.Aerodynamics.wing_loading.(caseName).dynamic_pressure = 0.5 * speeds(i)^2 * density;
end
end

