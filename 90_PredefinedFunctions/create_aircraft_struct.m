% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                             create_aircraft_struct.m
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
%       P R A K T I K U M   F L U G Z E U G E N T W U R F   1 9 / 20
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
% Creates the empty aircraft struct to be used for data exchange
% -------------------------------------------------------------------------
% DATE          VERSION     AUTHOR              CHANGELOG
% 15/10/2019    v1.0        Thomas Seren        initial release
% -------------------------------------------------------------------------
% Input variables:
% NAME          UNIT    DESCRIPTION
% aerodynamic_surface_cell - cell with the names of all aerodynamic surfaces
% variation_parameters - cell with the names of the variation parameters
% fix_parameters - cell with the names of the fixed parameters
% performance_information_fields - cell with the names of additional
%                                  performance parameters
% -------------------------------------------------------------------------
% Output variables:
% NAME          UNIT    DESCRIPTION
% Aircraft_Data -   struct with the aircraft definition
% -------------------------------------------------------------------------
% Limitations:  -
% Assumptions:  -
% -------------------------------------------------------------------------
% Example:      Aircraft_Data = create_aircraft_struct({'wing', 'HTP', 'VTP'}, {'wing_AR', 'wing_S_ref'}, {'TP_AR', 'static_margin'}, {'v_stall', 'v_cruise'})
% -------------------------------------------------------------------------
% Called functions:
% -
% -------------------------------------------------------------------------
% Error Management
% ERROR MESSAGE                            DESCRIPTION
% -
% -------------------------------------------------------------------------
% KNOWN ISSUES:
% - 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% first level struct description
function Aircraft_Data = create_aircraft_struct(aerodynamic_surface_cell, variation_parameters, fix_parameters, performance_information_fields)
    
    Aircraft_Data = struct(...
        'Mass', mass_struct(aerodynamic_surface_cell), ... 
        'Configuration', config_struct(variation_parameters, fix_parameters), ...
        'Performance', performance_struct(performance_information_fields), ...
        'Geometry', geometry_struct(aerodynamic_surface_cell), ...
        'Aerodynamics', aero_struct(), ...
        'Propulsion', propulsion_struct(), ...
        'Structure', structure_struct(), ...
        'id', NaN ...
        );
    printstructfields(Aircraft_Data) % create a text file in main work directory, that includes the complete struct overview
end
%% Mass Struct definition
function struct_out = mass_struct(aerodynamic_surface_cell)
    mass_table_variables = {'part','mass','cg_position','quantity'};
    
    struct_out = struct(...
        'total_mass', NaN, ... % [kg]
        'real_position_CG', NaN, ... % [m]
        'mass_table', array2table(zeros(20,numel(mass_table_variables)),'VariableNames', mass_table_variables) ...
        );
end
%% Configuration Struct definition
function struct_out = config_struct(variation_parameters, fix_parameters)
    struct_out = struct(...
        'design_parameter',cell2struct(cell(numel(variation_parameters),1), variation_parameters,1), ...
        'fixed_parameter', cell2struct(cell(numel(fix_parameters),1), fix_parameters,1) ...
        );
end
%% Perfomance Struct Definition
function struct_out = performance_struct(performance_information_fields)
    table_variables = {'flightphase', 'speed', 'time', 'current', 'power', 'energy','powersetting'};
    flightphases = {'take_off', 'climb', 'cruise straight','cruise turn', 'pattern straight', 'pattern turn'};
    performance_table = cell2table(cell(numel(flightphases), numel(table_variables)),'VariableNames', table_variables);
    performance_table(:,1) = flightphases.';
    for ii = 1:numel(performance_information_fields)
        struct_out.(performance_information_fields{ii}) = NaN;
    end
    struct_out.performance_table = performance_table;
end
%% Geometry Struct Definition
function struct_out = geometry_struct(aerodynamic_surface_cell)
    wing_table_variables = {'ID', 'x_offset', 'y_offset', 'z_offset', 'chord_length', 'twist'};
    for i_aero_surf = 1:numel(aerodynamic_surface_cell)
        struct_out.(aerodynamic_surface_cell{i_aero_surf})=struct(...
            'shape', cell2table(cell(1,numel(wing_table_variables)),'VariableNames', wing_table_variables), ...
            'airfoil', '', ...
            'projected_surface', NaN, ... % [m^2]
            'aspect_ratio', NaN, ... % [-]
            'taper_ratio', NaN, ... % [-]
            'mean_aerodynamic_chord', NaN, ... % [m]
            'position_neutral_point', NaN, ... % [m]
            'x_position', NaN, ... % [m]
            'volume', NaN, ... %[m^3]
            'wetted_surface', NaN, ... %[m^2]
            'thickness', NaN,... %[m]
            'span', NaN,... %[m]
            'root_aerodynamic_chord', NaN,... % [m]
            'dihedral',NaN... % [deg]
            );
    end
end
%% Aerodynamics Struct Definition
function struct_out = aero_struct()
    drag_table_variables = {'airspeed', 'glide_ratio', 'induced_drag', 'zero_lift_drag', 'total_drag','AoA'};
    table_variables = {'y_position', 'C_l', 'C_m', 'l_C_l'};
    wing_load_table = cell2table(cell(0,numel(table_variables)),'VariableNames', table_variables);
    relevant_speeds = {'stall', 'cruise', 'pattern', 'v_max_CL_max'};
    empty_wing_load_struct = struct(...
        'speed', NaN, ... % [m/s]
        'dynamic_pressure', NaN, ... % [Pa]
        'load_table', wing_load_table ...
        );
    for i_speed = 1:numel(relevant_speeds)
        wing_load_struct.(relevant_speeds{i_speed}) = empty_wing_load_struct;
    end
    struct_out = struct(...
        'target_position_CG', NaN, ... % [m]
        'position_neutral_point', NaN, ... % [m]
        'drag_coefficient_other', NaN, ... % [] accounts for drag from other components
        'speed_drag_table', cell2table(cell(6, numel(drag_table_variables)),'VariableNames', drag_table_variables), ...
        'wing_loading', wing_load_struct ...
        );
end
%% Propulsion Struct Definition
function struct_out = propulsion_struct()
    table_variables = {'speed', 'thrust', 'torque', 'rpm', 'throttle_setting', 'current', 'power', 'propeller_efficiency', 'motor_efficiency', 'overall_efficiency'};
    propulsion_table = cell2table(cell(1,numel(table_variables)),'VariableNames', table_variables);
    struct_out = struct(...
        'battery', struct('name', '', 'm', NaN, 'U_0', NaN, 'Capacity', NaN, 'Energy', NaN, 'C_Rate', NaN,'R_i', NaN, 'max_Current', NaN ), ...
        'motor', struct('name', '', 'm', NaN, 'I_idle', NaN, 'U_I_idle', NaN, 'R_i', NaN, 'kv', NaN, 'ke', NaN, 'kl', NaN,  'P_max', NaN), ...
        'propeller', struct('name', '', 'm', NaN, 'Diameter_m', NaN, 'Pitch_m', NaN, ...
            'dynamic_table', cell2table(cell(0,5),'VariableNames',{'J', 'CT', 'CP', 'CQ', 'eta'}), ...
            'static_table', cell2table(cell(0,4),'VariableNames',{'RPM', 'CT', 'CP', 'CQ'}) ...
            ), ...
        'full_throttle', propulsion_table, ...
        'part_throttle_18', propulsion_table, ...
        'part_throttle_19', propulsion_table, ...
        'part_throttle_20', propulsion_table, ...
        'part_throttle_21', propulsion_table, ...
        'part_throttle_22', propulsion_table ...
        );
        
end
%% Structure Struct Definition
function struct_out = structure_struct()
    table_variables = {'y_position', 'q', 'Q', 'M_b', 'M_t'};
    wing_forces_table = cell2table(cell(20,numel(table_variables)),'VariableNames', table_variables);    
    struct_out = struct(...
        'spar_diameter', NaN, ... % [m]
        'wing_twist', NaN, ... % [deg]
        'spar_web_thickness', NaN, ... % [m]
        'maximum_wing_load_table', wing_forces_table ...
        );
end
%% Function to create a text file in main work directory, that includes the complete struct overview
function printstructfields(struct_in)
    fid = fopen('aircraft_data_struct.txt','w');
    fprintf(fid,[inputname(1),'\n']);
    printstruct(struct_in,  1, fid);
    
    fclose(fid);
end
function printstruct(struct_in,n_tabs,fid)
    fieldnamecell = fieldnames(struct_in);
    for ii = 1:numel(fieldnamecell)
        fprintf(fid,['\n',repmat('\t',1,n_tabs),fieldnamecell{ii}]);
        if isstruct(struct_in.(fieldnamecell{ii}))
            printstruct(struct_in.(fieldnamecell{ii}), n_tabs + 1, fid);
        elseif istable(struct_in.(fieldnamecell{ii}))
            fprintf(fid,'   [ table: ');
            for variable = struct_in.(fieldnamecell{ii}).Properties.VariableNames
                fprintf(fid,variable{1});
                if ~isequal(variable,struct_in.(fieldnamecell{ii}).Properties.VariableNames(end))
                    fprintf(fid,' | ');
                end
            end
            fprintf(fid,' ]');
        elseif ~isnumeric(struct_in.(fieldnamecell{ii}))
            fprintf(fid,['   [ ',class(struct_in.(fieldnamecell{ii})),' ]']);
        end
        if n_tabs == 1
            fprintf(fid,'\n');
        end
    end
end
