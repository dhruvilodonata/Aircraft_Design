% Script to recalculate all combinations where a error occured which is not
% tide to the configuration itself
%% Initialising
clear; % delete everything from the workspace
close all % close any open figures for improved performance
warning off backtrace % make the warnings more beautiful
clc % clear command window

%% directories to be included in the workspace
data_directories.polars = './80_Databases/airfoil_Polars';
data_directories.airfoil_geometry = './80_Databases/airfoil_geometry';
data_directories.propulsion = './80_Databases/Propulsion';

%% Access Data
result_filename = 'Final_Configuration';
data = load(fullfile('70_Results',[result_filename,'.mat']));
combination_table = data.best_combination;
% indices = unique([find(isnan(combination_table.errorID));find(combination_table.errorID == 999)]);
% indices = find(combination_table.errorID == 999);
%indices = find(combination_table.points > 80);
indices = 1;
%combination_table = combination_table(indices,:);

%% Prepare Data
aerodynamic_surfaces = {'wing', 'tail'}; % define the names of all aerodynamic surfaces here
variation_parameter_names = combination_table.Properties.VariableNames(1:6);
fixed_parameter_names = combination_table.Properties.VariableNames(7:13);
performance_information_fields = combination_table.Properties.VariableNames(24:32);
airfoil_polar_naming_convention = '_T1_Re%f_M0.00_N5.0_XtrTop25%_XtrBot25%.txt'; % Naming convention coming from XFOIL/XFLR5
empty_aircraft_struct = create_aircraft_struct(aerodynamic_surfaces, variation_parameter_names,...
    fixed_parameter_names, performance_information_fields);

for i = 1:length(fixed_parameter_names)
    fixed_parameter.(fixed_parameter_names{i}) = unique(combination_table.(fixed_parameter_names{i})).';
end
for i = 1:length(variation_parameter_names)
    variation_parameter.(variation_parameter_names{i}) = unique(combination_table.(variation_parameter_names{i})).';
end
data_directories = rework_paths(data_directories);

%% Read Input Data
inputdata = read_inputdata(fixed_parameter,variation_parameter,data_directories,...
    airfoil_polar_naming_convention);

n_combinations = length(indices);
time_set_array = zeros(ceil(n_combinations/1000),1);
set_start = 1;
%% Calculate combinations
for i_set = set_start:ceil(n_combinations/1000)
    if i_set == ceil(n_combinations/1000)
        set = indices(((i_set-1)*1000+1):n_combinations);
    else
        set = indices(((i_set-1)*1000+1):((i_set)*1000));
    end
    text_msg(sprintf('Starting Configuration %d', (i_set-1)*1000+1));
    set_combination_table = combination_table(set,:);
    timer_set = tic;
    for i = 1:length(set)
        timer_iteration = tic;
        current_configuration = set_combination_table(i,:);
        aircraftdata = empty_aircraft_struct;
        for i_par = 1:numel(variation_parameter_names) % fill the aircraft struct with the configuration data
            aircraftdata.Configuration.design_parameter.(variation_parameter_names{i_par}) = current_configuration.(variation_parameter_names{i_par});
            if iscell(aircraftdata.Configuration.design_parameter.(variation_parameter_names{i_par}))
                aircraftdata.Configuration.design_parameter.(variation_parameter_names{i_par}) = aircraftdata.Configuration.design_parameter.(variation_parameter_names{i_par}){1};
            end
        end
        aircraftdata.id = current_configuration.id;
        for i_par = 1:numel(fixed_parameter_names) % fill the aircraft struct with the configuration data
            aircraftdata.Configuration.fixed_parameter.(fixed_parameter_names{i_par}) = current_configuration.(fixed_parameter_names{i_par});
            if iscell(aircraftdata.Configuration.fixed_parameter.(fixed_parameter_names{i_par}))
                aircraftdata.Configuration.fixed_parameter.(fixed_parameter_names{i_par}) = aircraftdata.Configuration.fixed_parameter.(fixed_parameter_names{i_par}){1};
            end
        end
        error_cell = cell(1,1);
        offsetdata = 0;
        errorID = NaN;
        try
            [aircraftdata,error_cell] = objective_function(aircraftdata,inputdata,error_cell);
        catch ME
            switch ME.identifier
                case 'Error999:FileID'
                    errorID = 999;
                case 'Error101:TrimError'
                    errorID = 101;
                case 'Error102:ReNumber'
                    errorID = 102;
                case 'Error103:CLmax'
                    errorID = 103;
                case 'Error104:stallSpeed'
                    errorID = 104;
                case 'Error201:maxThrust'
                    errorID = 201;
                case 'Error202:maxPower'
                    errorID = 202;
                case 'Error203:excessESCcurrent'
                    errorID = 203;
                case 'Error204:excessBatteryDiscargeCurrent'
                    errorID = 204;
                case 'Error301:noCG'
                    errorID = 301;
                case 'Error401:maxSpeed'
                    errorID = 401;
                case 'Error402:TOdistance'
                    errorID = 402;
                otherwise
                    warning(num2str(aircraftdata.id))
                    rethrow(ME)
            end
        end
        time_iteration = toc(timer_iteration);
        current_configuration.take_off_mass = aircraftdata.Mass.total_mass;
        if isnan(errorID)
            current_configuration.drag = aircraftdata.Aerodynamics.speed_drag_table.total_drag(5);
            current_configuration.tail_area = aircraftdata.Geometry.tail.projected_surface;
            current_configuration.tail_angle = aircraftdata.Geometry.tail.dihedral;
            current_configuration.tail_leverarm = aircraftdata.Geometry.tail.x_position;
            current_configuration.cg_position = aircraftdata.Mass.real_position_CG;
            current_configuration.position_neutral_point = aircraftdata.Aerodynamics.position_neutral_point;
            current_configuration.cruise_speed = aircraftdata.Performance.cruise_speed;
            current_configuration.climb_speed = aircraftdata.Performance.climb_speed;
            current_configuration.spar_web_thickness = aircraftdata.Structure.spar_web_thickness;
        end
        current_configuration.errorID = errorID;
        current_configuration.calculationTime = time_iteration;
    
        for ii = 1:numel(performance_information_fields)
            current_configuration.(performance_information_fields{ii}) = aircraftdata.Performance.(performance_information_fields{ii});
        end
        set_combination_table(i,:) = current_configuration;
    end
    time_set = toc(timer_set);
    if i_set < ceil(n_combinations/1000)
        time_set = toc(timer_set);
        time_set_array(i_set) = time_set;
        time_left = sum(time_set_array) / (i_set*1000) * (n_combinations - i_set*1000);
        if time_left > 24*3600
            time_finished = datestr(addtodate(now,time_left,'second'),'dd.mm, HH:MM');
        else
            time_finished = datestr(addtodate(now,time_left,'second'),'HH:MM');
        end
        text_msg(['The last set took ',num2str(round(time_set/60,0)),' minutes to calculate. Estimated time when finished is ',time_finished]);
    else
        text_msg(['The last set took ',num2str(round(time_set/60,0)),' minutes to calculate.']);
    end
    combination_table(set,:) = set_combination_table;
    save(fullfile('70_Results','combination_table_recalculation_recovery.mat'),'combination_table');
end
%save(fullfile('70_Results',[result_filename,'.mat']),'combination_table');