% Calculates the maximum possible speed in straight and level cruise flight
function [aircraftdata,error_cell] = maxSpeed(aircraftdata,inputdata,error_cell)
%% Extract values from input
S_ref = aircraftdata.Geometry.wing.projected_surface;
mass = aircraftdata.Mass.total_mass;

density = inputdata.density;
gravity = inputdata.gravity;

%% Calculation
for speed = 18:40
    q = 0.5 * density * speed^2;
    C_L = mass * gravity / (q * S_ref);
    [C_D, ~] = dragCoefficient(aircraftdata,inputdata,speed,C_L);
    drag = q * S_ref * C_D;
    [thrust,~,~,~] = powerCalculation(aircraftdata,inputdata,0,speed);
    if thrust < drag
        airspeed_max = speed;
        break
    elseif speed == 40
        airspeed_max = 41;
        warning('Maximum Airspeed higher than 40')
    end
end
if airspeed_max == 18
    error('Error401:maxSpeed','Maximum Speed is not high enough')
end

% For higher Fidelity:
if airspeed_max < 31
    for i = 1:10
        speed = airspeed_max - 0.1 * i;
        q = 0.5 * density * speed^2;
        C_L = mass * gravity / (q * S_ref);
        [C_D, ~] = dragCoefficient(aircraftdata,inputdata,speed,C_L);
        drag = q * S_ref * C_D;
        [thrust,~,~,~] = powerCalculation(aircraftdata,inputdata,0,speed);
        if thrust > drag
            airspeed_max = speed;
            break;
        end
    end
end
%% Write results to aircraft data struct
aircraftdata.Performance.maximum_speed = airspeed_max;
end