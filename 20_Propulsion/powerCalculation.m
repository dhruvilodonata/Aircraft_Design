% Calculates Power Consumption, current and throttle setting for given
% airspeed and thrust requirement. If Input Thrust is Zero, it ouputs the
% maximum thrust and power for given airspeed
function [thrust,power,current,throttle_setting] = powerCalculation(aircraftdata,inputdata,thrust,airspeed)

if any(strcmp(fieldnames(inputdata),'offsetdata'))
    offsetdata = inputdata.offsetdata;
    if any(strcmp(fieldnames(offsetdata),'power'))
        offsetFactor = offsetdata.power/100 + 1;
    else
        offsetFactor = 1;
    end
else
    offsetFactor = 1;
end

if length(thrust) ~= length(airspeed)
    if length(airspeed) == 1
        airspeed = ones(1,length(thrust)) * airspeed;
    elseif length(thrust) == 1
        thrust = ones(1,length(airspeed)) * thrust;
    else
        error('Length of thrust and airspeed do not match')
    end
end
full_throttle_table = aircraftdata.Propulsion.full_throttle;
i = 0;
part_throttle_table_cell = cell(1,5);
for v = 12:20
    i = i + 1;
    fieldName = ['part_throttle_',num2str(v)];
    part_throttle_table_cell{i} = aircraftdata.Propulsion.(fieldName);
end
power = zeros(1,length(thrust));
current = zeros(1,length(thrust));
throttle_setting = zeros(1,length(thrust));

for i = 1:length(thrust)
    if thrust(i) == 0
        thrust(i) = interp1(full_throttle_table.speed,full_throttle_table.thrust,airspeed(i));
        power(i) = offsetFactor * interp1(full_throttle_table.speed,full_throttle_table.power,airspeed(i));
        current(i) = interp1(full_throttle_table.speed,full_throttle_table.current,airspeed(i));
    else
        part_throttle_table = part_throttle_table_cell{airspeed(i) - 11};
        if thrust > max(part_throttle_table.thrust)
            power(i) = NaN;
            current(i) = NaN;
            throttle_setting(i) = NaN;
        else
            power(i) = offsetFactor * interp1(part_throttle_table.thrust,part_throttle_table.power,thrust(i));
            current(i) = interp1(part_throttle_table.thrust,part_throttle_table.current,thrust(i));
            throttle_setting(i) = interp1(part_throttle_table.thrust,part_throttle_table.throttle_setting,thrust(i));
        end
    end
end
end