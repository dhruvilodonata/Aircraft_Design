clear
clc
load(fullfile('70_Results','Final_Configuration.mat'))
error_cell = cell(1,1);
motor_ids = 1:5;
propeller_ids = 1:4;
battery_ids = 1:3;
n_combinations = length(motor_ids)*length(propeller_ids)*length(battery_ids);
T = array2table(zeros(n_combinations,6),'VariableNames',{'points','propeller','motor','battery','dragCruise','dragPattern'});
n = 0;
%% Calculation
for i = 4%1:length(propeller_ids)
    for j = 2%1:length(motor_ids)
        for k = 3%1:length(battery_ids)
            n = n + 1;
            aircraftdata_new = aircraftdata;
            aircraftdata_new.Performance.stall_speed = NaN;
            aircraftdata_new.Configuration.design_parameter.propeller_ID = propeller_ids(i);
            aircraftdata_new.Configuration.design_parameter.motor_ID = motor_ids(j);
            aircraftdata_new.Configuration.design_parameter.battery_ID = battery_ids(k);
            
            try
                [aircraftdata_new,error_cell] = objective_function(aircraftdata_new,inputdata,error_cell);
                T.points(n) = aircraftdata_new.Performance.points;
                T.dragPattern(n) = aircraftdata_new.Aerodynamics.speed_drag_table.total_drag(5);
                T.dragCruise(n) = aircraftdata_new.Aerodynamics.speed_drag_table.total_drag(3);
                T.propeller(n) = propeller_ids(i);
                T.motor(n) = motor_ids(j);
                T.battery(n) = battery_ids(k);
            catch
            end
        end
    end
end
