function [ aircraftdata, error_cell ] = Performance( aircraftdata, inputdata, error_cell)
% function objectives:
% 1) calculate flight time (from which pattern_time is needed for points function)
% 2) output performance table
% 3) stall speed and take-off distance

%% Define Fixed Parameters
radius_cruise = 50; % [m] radius for cruise laps
radius_pattern = 40; % [m] radius for reconnaissance pattern
power_avionics = 7; % [W] Power Usage of Avionics
cruise_altitude_AGL = 100; % [m] Cruise Altitude Above Ground Level
friction_coeff = 0.08; % rolling friction co efficient of wheel on ground

%% Extract Input Data
g = inputdata.gravity;
density = inputdata.density;
speed_pattern = inputdata.speed_pattern;

%% Extract from Aircraft Data Struct
stall_speed = aircraftdata.Performance.stall_speed; % !written by Aero team
mass = aircraftdata.Mass.total_mass; % !written by Structure team
S_ref = aircraftdata.Configuration.design_parameter.wing_S_ref; % !written by Aero team
energy_MAX = 3600e-3 * aircraftdata.Propulsion.battery.Q * aircraftdata.Propulsion.battery.U_nom * 0.8;    % in J

%% Take off simulation
C_L_TO = CLforMinCD(aircraftdata,inputdata,8);
time = 0;
distance = 0;
velocity = 0;
acceleration = 0;
timeIncrement = 0.1;
energy_TO = 0;
airspeed_TO = stall_speed + 3;
AoA = 5;    % AoA in degree when rolling on the ground
while true
    distance = distance + velocity * timeIncrement;
    velocity = velocity + acceleration * timeIncrement;
    [thrust_TO,power_TO,current_TO,powersetting_TO] = ...
        powerCalculation(aircraftdata,inputdata,0,velocity);
    energy_TO = energy_TO + (power_TO + power_avionics) * timeIncrement;
    if velocity < 3
        [C_D_TO,~] = dragCoefficient(aircraftdata,inputdata,velocity,0);
    else
        [C_D_TO,~] = dragCoefficient(aircraftdata,inputdata,velocity,C_L_TO);
    end
    lift_TO = 0.5 * S_ref*density*velocity^2*C_L_TO;
    drag_aero = 0.5 * S_ref*density*velocity^2*C_D_TO;
    drag_friction = friction_coeff*(mass*g-lift_TO);
    drag_TO = drag_aero + drag_friction;
    acceleration = (thrust_TO*cosd(AoA)-drag_TO)/mass;
    time = time + timeIncrement;
    if velocity >= airspeed_TO
        break
    elseif distance > 50
        error('Error402:TOdistance','No T/O inside 50m')
    end
end
speed_TO = velocity;
time_TO = time;
distance_TO = distance;
C_L_TO = 2 * mass * g * 1.2 / (density * S_ref * speed_TO^2);
[C_D_TO,AoA_TO] = dragCoefficient(aircraftdata,inputdata,speed_TO, C_L_TO);
drag_TO = 0.5*S_ref*density*velocity^2*C_D_TO;
glideRatio_TO = C_L_TO / C_D_TO;

%% Climb rate and time to reach 100 m height with power required 
% checking best AoC = Angle of Climb and climb speed combination to reach in min time 
AoC = 5 * pi / 180; % [radian] start from 1 degree
velocities_climb = 12:20;
energy_climb_array = zeros(1,length(velocities_climb))*inf;
time_climb_array = ones(1,length(velocities_climb))*inf;
power_climb_array = zeros(1,length(velocities_climb));
drag_climb_array = zeros(1,length(velocities_climb));
glideRatio_climb_array = zeros(1,length(velocities_climb));
current_climb_array = zeros(1,length(velocities_climb));
powersetting_climb_array = zeros(1,length(velocities_climb));
AoA_climb = zeros(1,length(velocities_climb));
for i = 1:length(velocities_climb)
    speed_climb = velocities_climb(i);
    difference = 100;
    j = 0;
    while difference > 5 && AoC > 0
        j = j + 1;
        C_L_climb = 2*mass*g*cos(AoC)/(speed_climb^2*S_ref*density);
        [C_D_climb,AoA_climb(i)] = dragCoefficient(aircraftdata,inputdata,speed_climb,C_L_climb);
        [thrust_climb,power_climb_array(i),current_climb_array(i),powersetting_climb_array(i)] =...
            powerCalculation(aircraftdata,inputdata,0,speed_climb);
        drag_climb_array(i) = (density*S_ref*speed_climb^2*C_D_climb)/2;
        glideRatio_climb_array(i) = C_L_climb / C_D_climb;
        epsilon = C_D_climb/C_L_climb; %glide angle
        AoC_new = asin(thrust_climb/(mass*g)/sqrt(1+epsilon^2)) - atan(epsilon); % [radian] AoA has to be taken into account
        %AoC_2 = asin((thrust-drag_climb)/mass*g); % [radian]
        difference = abs((AoC - AoC_new)*100/AoC);
        AoC = AoC_new;
        if j > 5
            power_climb_array(i) = inf;
            AoC = 0;
            break
        end
    end
    RoC = speed_climb*sin(AoC); % rate of climb
    if RoC > 0
        time_climb_array(i) = cruise_altitude_AGL/RoC; %[sec]
    else
        time_climb_array(i) = inf;
    end
    energy_climb_array(i) = (power_climb_array(i) + power_avionics)*time_climb_array(i);
end
if ~any(energy_climb_array ~= inf)
    error('Error403:climbPower','Not enough Power to achieve positive Climb rate')
end
[energy_climb,index] = min(energy_climb_array);
powersetting_climb = powersetting_climb_array(index);
speed_climb = velocities_climb(index);
time_climb = time_climb_array(index);
power_climb = power_climb_array(index);
drag_climb = drag_climb_array(index);
glideRatio_climb = glideRatio_climb_array(index);
current_climb = current_climb_array(index);
AoA_climb = AoA_climb(index);

%% Cruise Calculation
speed_cruise = 13:18;
energy_cruise = zeros(length(speed_cruise),2) * inf;
power_cruise = zeros(length(speed_cruise),2);
drag_cruise = zeros(length(speed_cruise),2);
glideRatio_cruise = zeros(length(speed_cruise),2);
current_cruise = zeros(length(speed_cruise),2);
time_cruise = zeros(length(speed_cruise),2);
powersetting_cruise = zeros(length(speed_cruise),2);
AoA_cruise = zeros(length(speed_cruise),2);
for i = 1:length(speed_cruise)
    distance_straight_cruise_lap = 2 * 200; %[m]
    time_straight_cruise = 10 * distance_straight_cruise_lap / speed_cruise(i);
    distance_turn_cruise = 2 * pi * radius_cruise;
    time_turn_cruise = 10 * distance_turn_cruise / speed_cruise(i);
    time_cruise(i,:) = [time_straight_cruise, time_turn_cruise];

    loadFactor_cruise_straight = 1;
    loadFactor_cruise_turn = sqrt((speed_cruise(i)^4/(radius_cruise * g)^2)+1);
    loadFactor_cruise = [loadFactor_cruise_straight, loadFactor_cruise_turn];

    C_L_cruise = 2 * mass * g * loadFactor_cruise / (density * S_ref * speed_cruise(i)^2);
    [C_D_cruise,AoA_cruise(i,:)] = dragCoefficient(aircraftdata,inputdata,speed_cruise(i),C_L_cruise);
    drag_cruise(i,:) = 0.5 * S_ref * density * speed_cruise(i)^2 * C_D_cruise;
    glideRatio_cruise(i,:) = C_L_cruise ./ C_D_cruise;
    thrust_cruise = drag_cruise(i,:) ./ cosd(AoA_cruise(i,:));
    [~,power_cruise(i,:),current_cruise(i,:),powersetting_cruise(i,:)] =...
        powerCalculation(aircraftdata,inputdata,thrust_cruise,speed_cruise(i));
    if any(isnan(sum(power_cruise,2)))
        energy_cruise(i,:) = inf;
        break
    end
    energy_cruise(i,:) = power_cruise(i,:) .* time_cruise(i,:) + power_avionics * time_cruise(i,:);
end
if ~any(sum(energy_cruise,2) ~= inf)
    error('Error404:cruisePower','Not enough Power to sustain Cruise Flight')
end
[~, index] = min(sum(energy_cruise,2));
time_cruise = time_cruise(index,:);
energy_cruise = energy_cruise(index,:);
speed_cruise = speed_cruise(index);
power_cruise = power_cruise(index,:);
drag_cruise = drag_cruise(index,:);
glideRatio_cruise = glideRatio_cruise(index,:);
current_cruise = current_cruise(index,:);
powersetting_cruise = powersetting_cruise(index,:);
AoA_cruise = AoA_cruise(index,:);

%% Pattern Calculation
% 1 Lap Pattern
distance_straight_pattern = 6 * 300 + 320; % [m]
time_straight_pattern = distance_straight_pattern / speed_pattern;
distance_turn_pattern = 3 * 2 * pi * radius_pattern;
time_turn_pattern = distance_turn_pattern / speed_pattern;
time_lap_pattern = [time_straight_pattern, time_turn_pattern];

% Energy left for Pattern flight
energy_used = energy_TO + energy_climb + sum(energy_cruise);
energy_pattern_overall = energy_MAX - energy_used;

loadFactor_pattern_straight = 1;
loadFactor_pattern_turn = sqrt((speed_pattern^4/(radius_pattern * g)^2)+1);
loadFactor_pattern = [loadFactor_pattern_straight, loadFactor_pattern_turn];

C_L_pattern = 2 * mass * g * loadFactor_pattern / (density * S_ref * speed_pattern^2);
[C_D_pattern,AoA_pattern] = dragCoefficient(aircraftdata,inputdata,speed_pattern,C_L_pattern);
glideRatio_pattern = C_L_pattern ./ C_D_pattern;
drag_pattern = 0.5 * S_ref * density * speed_pattern^2 * C_D_pattern;
thrust_pattern = drag_pattern ./ cosd(AoA_pattern);
[~,power_pattern,current_pattern,powersetting_pattern] =...
    powerCalculation(aircraftdata,inputdata,thrust_pattern,speed_pattern);
if any(isnan(power_pattern))
    error('Error201:maxThrust','Not enough thrust to sustain pattern flight')
end
energy_lap_pattern = power_pattern .* time_lap_pattern + power_avionics * time_lap_pattern;

number_laps_pattern = energy_pattern_overall / sum(energy_lap_pattern);
time_pattern = time_lap_pattern * number_laps_pattern;
energy_pattern = energy_lap_pattern * number_laps_pattern;

time_overall = time_TO + time_climb + sum(time_cruise) + sum(time_pattern);

%% Write Output to Aircraft Data Struct

aircraftdata.Performance.performance_table.speed = [speed_TO, speed_climb,speed_cruise,speed_cruise,speed_pattern,speed_pattern]';
aircraftdata.Performance.performance_table.time = [time_TO,time_climb,time_cruise,time_pattern]'; %?
aircraftdata.Performance.performance_table.current = [current_TO,current_climb,current_cruise,current_pattern]'; %?
aircraftdata.Performance.performance_table.power = [power_TO,power_climb,power_cruise,power_pattern]';
aircraftdata.Performance.performance_table.energy = [energy_TO,energy_climb,energy_cruise,energy_pattern]';
aircraftdata.Performance.performance_table.powersetting = [powersetting_TO,powersetting_climb,powersetting_cruise,powersetting_pattern]';

aircraftdata.Performance.take_off_distance = distance_TO;
aircraftdata.Performance.pattern_time = sum(time_pattern);
aircraftdata.Performance.overall_flighttime = time_overall;
aircraftdata.Performance.cruise_speed = speed_cruise;
aircraftdata.Performance.climb_speed = speed_climb;

aircraftdata.Aerodynamics.speed_drag_table.airspeed = [speed_TO, speed_climb, speed_cruise, speed_cruise, speed_pattern, speed_pattern]';
aircraftdata.Aerodynamics.speed_drag_table.glide_ratio = [glideRatio_TO, glideRatio_climb, glideRatio_cruise, glideRatio_pattern]';
aircraftdata.Aerodynamics.speed_drag_table.total_drag = [drag_TO, drag_climb, drag_cruise, drag_pattern]';
aircraftdata.Aerodynamics.speed_drag_table.AoA = [AoA_TO, AoA_climb, AoA_cruise, AoA_pattern]';
end