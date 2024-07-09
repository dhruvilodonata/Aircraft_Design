function [ aircraftdata, error_cell ] = Points( aircraftdata, inputdata, error_cell)

time_pattern = aircraftdata.Performance.pattern_time;
%% fixed setups
print_time_max = 110*60*60 ;

%% stall speed and take-off distance penalty factor a
a = 1;

if aircraftdata.Performance.stall_speed > 9
    %a = a * 0.5;
end

if aircraftdata.Performance.take_off_distance > 40
    a = a * 0.8;
end

%% P for points
if time_pattern <= 0
    P_flyoff = 0;
else
    P_flyoff = 12/60 * time_pattern;
end
%P_prediction = 16; % set to 80% points temporarily (10% error in predicting flyoff time)
P_print = 100 * (1 - aircraftdata.Performance.print_time / print_time_max ); % unit [s]
% reminder: put error warning when print_time exceeds max @ STRUCTURE team

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_total = a * P_flyoff + P_print;
aircraftdata.Performance.points = P_total;
end

