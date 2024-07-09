% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                             Run.m
%
%      L E H R S T U H L   F U E R   L U F T F A H R T S Y S T E M E
%        I N S T I T U T E   O F   A I R C R A F T   D E S I G N
%                      Web: www.lls.lrg.tum.de
%                        ______
%                          |  |    |\  /|
%                          |  |    | \/ |
%                          |  |____|    |
%                Technische Universitaet Muenchen TUM
%
% (c) 2019 by Institute of Aircraft Design of TUM
%                       All Rights Reserved
%
%       P R A K T I K U M   F L U G Z E U G E N T W U R F 
%                             __   __   __   __ 
%                        |   |__| |__| |__| |__|
%                        |   |  | |\   |    |  |
%                        |__ |  | | \  |    |  |
%                        
%                 L L S   Advanced Research Project Agency
%
% Enter:        -
% Module:       -
% Version:      1.0
% -------------------------------------------------------------------------
% DESCRIPTION:
% Matlab Script to run a design parameter variation and mass iteration for
% the pracitcal course Aircraft Design.
% -------------------------------------------------------------------------
% DATE          VERSION     AUTHOR              CHANGELOG
% 28/10/2019    v1.0        Thomas Seren        initial release
% -------------------------------------------------------------------------
% Limitations:  -
% Assumptions:  -
% -------------------------------------------------------------------------
% Example:      Run.m
% -------------------------------------------------------------------------
% Called functions:
% -
% -------------------------------------------------------------------------
% Error Management
% ERROR MESSAGE                            DESCRIPTION
% -                                         -
% -------------------------------------------------------------------------
% KNOWN ISSUES:
% - 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Initialising
clear; % delete everything from the workspace
close all % close any open figures for improved performance
warning off backtrace % make the warnings more beautiful
clc % clear command window
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% modify to your needs %%%%%%%%%%%%%%%%%%%%%%%%%

%% start timer
globalTimer = tic;

%% directories to be included in the workspace

data_directories.polars = './80_Databases/airfoil_Polars';
data_directories.airfoil_geometry = './80_Databases/airfoil_geometry';
data_directories.propulsion = './80_Databases/Propulsion';

sub_directories = {...
    './10_Aerodynamics'
    './20_Propulsion'
    './30_Structure'
    './40_Performance'
    './50_OtherFunctions'
    './70_Results'
    './80_Databases'
    './90_PredefinedFunctions'
    './99_MatlabFileExchange'
    ... % add all directories that you need
    }; 

%% Design parameter ranges:

% Every design/variation parameter that is worth to be looked into needs to
% get an vector (for numeric values) or cell (for text values) assigned.
% update this to your needs
            
variation_parameter.wing_S_ref = 0.225:0.001:0.235; % [m^2]
variation_parameter.wing_aspectRatio = 12:1:17; % [-]
variation_parameter.wing_airfoil = {'CLARK_Y','FX_60_126','FX_63_137','NACA_6412'};

variation_parameter.motor_ID = 1:5;
variation_parameter.propeller_ID = 1:4;
variation_parameter.battery_ID = 2:3;
% ...

%% Invalid Motor - Propeller Combinations

% Not all combinations of propeller and motors are feasible, due to the
% maximum power a motor is allowed to use - by selecting which combinations
% are invalid, the simulation time can be drastically reduced!
% mx2 array - first number motor_ID, second number propeller_ID

invalid_motor_propeller_combinations = [...
    ];


%% fixed parameter definition

aerodynamic_surfaces = {'wing', 'tail'}; % define the names of all aerodynamic surfaces here
airfoil_polar_naming_convention = '_T1_Re%f_M0.00_N5.0_XtrTop25%_XtrBot25%.txt'; % Naming convention coming from XFOIL/XFLR5

% insert all parameters here, that need to be fixed, to completely define
% the aircraft. The following are only examples and not close to being
% complete
fixed_parameter.spar_position = 0.25; % [-]
fixed_parameter.wing_taperRatio = 0.5;
fixed_parameter.tail_aspectRatio = 4;
fixed_parameter.tail_taperRatio = 0.6;
fixed_parameter.dihedral = 2;
fixed_parameter.dihedral_position = 0;
fixed_parameter.tail_airfoil = 'NACA_0012';

%%%%%%%%%%%%%%%%%%%%%%%%%%% modify till here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% information needed for a quick evaluation of the configuration
additional_fields = {'calculationTime', 'errorID', 'take_off_mass','id','drag','tail_area','tail_angle','tail_leverarm','cg_position','position_neutral_point','spar_web_thickness'};
performance_information_fields = {'points', 'overall_flighttime', 'take_off_distance', 'stall_speed', 'cruise_speed', 'pattern_time', 'print_time', 'maximum_speed', 'climb_speed'};

%% add paths to working directory

for ii = 1:numel(sub_directories)
    addpath(genpath(sub_directories{ii}));
end
%% Starting message

dispheader()
text_msg('Start');

%% Definition of some variables.

variation_parameter_names = fieldnames(variation_parameter).';
fixed_parameter_names = fieldnames(fixed_parameter).';
data_directories = rework_paths(data_directories);
sub_directories = rework_paths(sub_directories);

datestring = datestr(now,'YYYYmmddHHMM'); % to be used to store data from different runs

error_cell = cell(0,1); % to be used to store all errors/warnings that might happen

%% reshape all variation vectors, in a 1xn dimension

for i_par = 1:numel(variation_parameter_names)
    variation_parameter.(variation_parameter_names{i_par}) = reshape(variation_parameter.(variation_parameter_names{i_par}),1,[]);
end

%%
text_msg('Read Inputs');
%% Read input data
inputdata = read_inputdata(fixed_parameter,variation_parameter,data_directories,...
    airfoil_polar_naming_convention);

%% Check if any errors occured during the preparation and stop the process
if numel(error_cell) > 0
    text_msg('Errors Detected')
    error_msg('ENDING THE SCRIPT, BECAUSE THERE ARE ERRORS!');
%     cleanworkspace(sub_directories)
%     return % end the script if any errors/warnings have occured
end

% continue if there haven't been any errors
text_msg('Inputs Successfully Read'); 
%% Creation Design Parameter Combination Table
for i_par = 1:numel(variation_parameter_names)
    values = variation_parameter.(variation_parameter_names{i_par});
    if iscell(values)
        values = 1:numel(values); % can't combine text values
    end
    if i_par == 1
        combination_array = values; % first vector is the initial vector
    else
        combination_array = combvec(combination_array,values); % create combination of current array with new vector
    end
end
combination_array = combination_array.'; % transpose the array

% invalid prop-motor combination deselection
deselection_vector = false(size(combination_array,1),1); % create an initial deselction array for sorting out the invalid motor-propeller combinations
for ii = 1:size(invalid_motor_propeller_combinations,1)
    deselection_vector = deselection_vector | ...
        (...
        combination_array(:,cellcmp(variation_parameter_names,'motor_ID')) == invalid_motor_propeller_combinations(ii, 1) ... % find all motors with given ID
        & ...
        combination_array(:,cellcmp(variation_parameter_names,'propeller_ID')) == invalid_motor_propeller_combinations(ii, 2) ... % find all propellers with given ID
        );
end
combination_array(deselection_vector,:) = []; % clear all rows with invalid motor-propeller combinations

if isfile(fullfile('70_Results','combination_table_recovery.mat'))
    while true
        recover_str = input("Do you want to use the recovered Combination Table? (y/n)",'s');
        if strcmp(recover_str,'yes') || strcmp(recover_str,'y')
            recover = true;
            break
        elseif strcmp(recover_str,'no') || strcmp(recover_str,'n')
            recover = false;
            break
        else
            warning('wrong user input')
        end
    end
else
    recover = false;
end
if recover
    recovered_data = load(fullfile('70_Results','combination_table_recovery'));
    combination_table = recovered_data.combination_table;
    set_start = ceil(find(isnan(combination_table.id), 1, 'first') / 1000);
else
    set_start = 1;
    % table creation
    combination_table = array2table(combination_array, 'VariableNames', variation_parameter_names); % turn Array into table for better readability
    
    % replace numbers with text, where applicable:
    for i_par = 1:numel(variation_parameter_names)
        if iscell(variation_parameter.(variation_parameter_names{i_par}))
            combination_table.(variation_parameter_names{i_par}) = (variation_parameter.(variation_parameter_names{i_par})(combination_table.(variation_parameter_names{i_par}))).'; 
        end
    end
    
     % add fixed parameter
    for i_par = 1:numel(fixed_parameter_names)
        fieldname = fixed_parameter_names{i_par};
        if ischar(fixed_parameter.(fieldname))
            combination_table.(fieldname) = repmat({fixed_parameter.(fieldname)},height(combination_table),1);
        else
            combination_table.(fieldname) = ones(height(combination_table),1) * fixed_parameter.(fieldname);
        end
    end
    % add additonal fields
    combination_table = add_NaN(combination_table, additional_fields);
    
    % add performance information fields
    combination_table = add_NaN(combination_table, performance_information_fields);
end

%% Parameter variation
n_combinations = height(combination_table); % number of combinations
time_set_array = zeros(ceil(n_combinations/1000),1);
%% Create empty Aircraft Struct
empty_aircraft_struct = create_aircraft_struct(aerodynamic_surfaces, variation_parameter_names, fixed_parameter_names, performance_information_fields);
%%
text_msg(sprintf('Calculating %d Configurations', n_combinations));
for i_set = set_start:ceil(n_combinations/1000)
    if i_set == ceil(n_combinations/1000)
        set = ((i_set-1)*1000+1):n_combinations;
    else
        set = ((i_set-1)*1000+1):((i_set)*1000);
    end
    text_msg(sprintf('Starting Configuration %d', (i_set-1)*1000+1));
    set_combination_table = combination_table(set,:);
    timer_set = tic;
    parfor i_combination = 1:length(set) % iterates through all configurations
        timer_iteration = tic;
        error_cell = cell(0,1); % to be used to store all errors/warnings that might happen
        current_configuration = set_combination_table(i_combination,:);
        aircraftdata = empty_aircraft_struct;
        %%
        for i_par = 1:numel(variation_parameter_names) % fill the aircraft struct with the configuration data
            aircraftdata.Configuration.design_parameter.(variation_parameter_names{i_par}) = current_configuration.(variation_parameter_names{i_par});
            if iscell(aircraftdata.Configuration.design_parameter.(variation_parameter_names{i_par}))
                aircraftdata.Configuration.design_parameter.(variation_parameter_names{i_par}) = aircraftdata.Configuration.design_parameter.(variation_parameter_names{i_par}){1};
            end
        end
        aircraftdata.id = i_combination + (i_set-1)*1000;
        for i_par = 1:numel(fixed_parameter_names) % fill the aircraft struct with the configuration data
            aircraftdata.Configuration.fixed_parameter.(fixed_parameter_names{i_par}) = current_configuration.(fixed_parameter_names{i_par});
            if iscell(aircraftdata.Configuration.fixed_parameter.(fixed_parameter_names{i_par}))
                aircraftdata.Configuration.fixed_parameter.(fixed_parameter_names{i_par}) = aircraftdata.Configuration.fixed_parameter.(fixed_parameter_names{i_par}){1};
            end
        end
        errorID = NaN;
        try
            offsetdata = 0;
            aircraftdata = objective_function(aircraftdata,inputdata,offsetdata);
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
                case 'Error204:excessBatteryDischargeCurrent'
                    errorID = 204;
                case 'Error301:noCG'
                    errorID = 301;
                case 'Error302:printingTime'
                    errorID = 302;
                case 'Error303:safetyFactor'
                    errorID = 303;
                case 'Error304:twistAngle'
                    errorID = 303;
                case 'Error401:maxSpeed'
                    errorID = 401;
                case 'Error402:TOdistance'
                    errorID = 402;
                case 'Error403:climbPower'
                    errorID = 403;
                case 'Error404:cruisePower'
                    errorID = 404;
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
        current_configuration.id = aircraftdata.id;
    
        for ii = 1:numel(performance_information_fields)
            current_configuration.(performance_information_fields{ii}) = aircraftdata.Performance.(performance_information_fields{ii});
        end
        set_combination_table(i_combination,:) = current_configuration;
    
    end % End Design Parameter-Loop
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
    save(fullfile('70_Results','combination_table_recovery.mat'),'combination_table');
end
%% save data
text_msg('Save Results')
save(fullfile('70_Results', [datestring,'_combination_table.mat']),'combination_table');
delete(fullfile('70_Results','combination_table_recovery.mat'));
%% End timer
time_iteration = toc(globalTimer);

%% End of script
text_msg('All Calculations done')
% cleanworkspace(sub_directories)

%% subfunctions
        
%% Clean Workspace:
% remove sub directories from working path
function cleanworkspace(sub_directories)
    warning on backtrace
    for ii = 1:numel(sub_directories)
        rmpath(genpath(sub_directories{ii}));
    end
    disp('****** End                                                                        ******');
    for ii=1:4
        disp('****************************************************************************************');
    end
end
%%
function ia = cellcmp(cell_array, string)
    [~, ia] = intersect(cell_array, string);
end
%%
function combination_table = add_NaN(combination_table, additional_fields)
    for field = additional_fields % add additional output parameters, that need to be written
        combination_table.(field{1}) = NaN(height(combination_table),1);
    end
end
%%
function dispheader()
    disp('****************************************************************************************');
    disp('****************************************************************************************');
    disp('****************************************************************************************');
    disp('****************************************************************************************');
    disp('******                                                                            ******');
    disp('******        I N S T I T U T E   O F   A I R C R A F T   D E S I G N             ******');
    disp('******                      Web: www.lls.lrg.tum.de                               ******');
    disp('******                        ______                                              ******');
    disp('******                          |  |    |\  /|                                    ******');
    disp('******                          |  |    | \/ |                                    ******');
    disp('******                          |  |____|    |                                    ******');
    disp('******                 Technichal University of Munich TUM                        ******');
    disp('******                                                                            ******');
    disp('******        (c) 2020 by Institute of Aircraft Design of TUM                     ******');
    disp('******                       All Rights Reserved                                  ******');
    disp('******                                                                            ******');
    disp('******      P R A C T I C A L   C O U R S E:   A I R C R A F T   D E S I G N      ******');
    disp('******                             __   __   __   __                              ******');
    disp('******                        |   |__| |__| |__| |__|                             ******');
    disp('******                        |   |  | |\   |    |  |                             ******');
    disp('******                        |__ |  | | \  |    |  |                             ******');
    disp('******                                                                            ******');
    disp('******               L L S   Advanced Research Project Agency                     ******');
    disp('******                                                                            ******');
    disp('****************************************************************************************');
    disp('****************************************************************************************');
    disp('****************************************************************************************');
    disp('****************************************************************************************');

end