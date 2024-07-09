%% Preparation
clear
clc
load(fullfile('70_Results','Final_Configuration.mat'));
%% Initialisation
airspeed = 8.9:0.1:35;
drag = zeros(1,length(airspeed));
ClCd = zeros(1,length(airspeed));
%% Reading of Inputs
mass = aircraftdata.Mass.total_mass;
S_ref = aircraftdata.Geometry.wing.projected_surface;
density = inputdata.density;
gravity = inputdata.gravity;
%% Calculation
for i = 1:length(airspeed)
    C_L = 2 * mass * gravity / (airspeed(i)^2 * S_ref * density);
    try
        [C_D, AoA] = dragCoefficient(aircraftdata,inputdata,airspeed(i),C_L);
    catch ME
        switch ME.identifier
            case 'Error105:dragStall'
                C_D = NaN;
            otherwise
                rethrow(ME)
        end
    end
    drag(i) = 0.5 * airspeed(i)^2 * S_ref * density * C_D;
    ClCd(i) = C_L / C_D;
end
T = array2table([airspeed',drag',ClCd'],'VariableNames',{'Airspeed','Drag','CL/CD'});