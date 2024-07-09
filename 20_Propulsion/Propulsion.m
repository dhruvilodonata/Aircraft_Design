% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                            Propulsion.m
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
% (c) 2022 by Institute of Aircraft Design of TUM
%                       All Rights Reserved
%
%       P R A K T I K U M   F L U G Z E U G E N T W U R F   21 / 22
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
% Main propulsion function that fills the aircraftdata struct with the data
% of the propulsion system at a given flight state (airspeed)
% -------------------------------------------------------------------------
% DATE          VERSION     AUTHOR              CHANGELOG
% 12/01/22    v1.0        Fabian Zach        initial release
% -------------------------------------------------------------------------
% Input variables:
% NAME                      UNIT        DESCRIPTION
% aircraftdata              [-]         aircraftdata struct
% Propulsion_Data_Struct    [-]         struct containing the propulsion
%                                       data
% propulsion_lookup         [-]         struct containig the lookup tables
%                                       for the propulsion configurations
% ID_P                      [-]         Propeller ID
% ID_M                      [-]         Motor ID
% ID_B                      [-]         Battery ID
% airspeed                  [m s^-1]    Expected airspeed 
% error_cell                [-]         Error cell
% -------------------------------------------------------------------------
% Output variables:
% NAME                  UNIT        DESCRIPTION
% aircraftdata          [-]         New aircraftdata struct
% error_cell            [-]         New error cell
% -------------------------------------------------------------------------
% Limitations:  -
% Assumptions:  -
% -------------------------------------------------------------------------
% Example:      [interp_lookup_table] =
% interpolate_lookup(Propulsion_Data_Struct,1,1,1,18)
% -------------------------------------------------------------------------
% Called functions:
% -
% -------------------------------------------------------------------------
% Error Management
% ERROR MESSAGE                            DESCRIPTION
% 
% -------------------------------------------------------------------------
% KNOWN ISSUES:
% - 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [ aircraftdata, error_cell ] = Propulsion(aircraftdata, inputdata, error_cell)
    % Required input data: thrust, component IDs
    ID_P = aircraftdata.Configuration.design_parameter.propeller_ID;
    ID_M = aircraftdata.Configuration.design_parameter.motor_ID;
    ID_B = aircraftdata.Configuration.design_parameter.battery_ID;
    %airspeed = inputdata.speed_pattern;
    Propulsion_Data_Struct = inputdata.propulsion.Propulsion_Data_Struct;
    propulsion_lookup = inputdata.propulsion.propulsion_lookup;
    
    i = 0;
    for airspeed = 12:20
        i = i + 1;
        % Interpolate complete lookup table at given airspeed
        interp_lookup_table = interpolate_lookup(propulsion_lookup,ID_P,ID_M,ID_B,airspeed);
        T = interp_lookup_table.T;
        Q = interp_lookup_table.Q;
        RPM = interp_lookup_table.RPM;
        U_motor = interp_lookup_table.V_motor;
        U_ESC = interp_lookup_table.V_ESC;
        I_ESC = interp_lookup_table.I_ESC;
        P_tot = interp_lookup_table.P_tot;
        I_motor = interp_lookup_table.I_motor;
        P_motor = interp_lookup_table.P_motor;
        eta_prop = interp_lookup_table.eta_prop;
        eta_motor = interp_lookup_table.eta_motor;
        eta_tot = interp_lookup_table.eta_tot;
        airspeed_vec = airspeed * ones(max(size(T)),1);
    
        % Generate the part_throttle table for the aircraftdata struct
        fieldName = ['part_throttle_',num2str(airspeed)];
        aircraftdata.Propulsion.(fieldName) = table(airspeed_vec,T,Q,RPM,U_motor./U_ESC, I_ESC, P_tot, I_motor, P_motor, eta_prop, eta_motor,eta_tot,'VariableNames',{'speed' 'thrust' 'torque' 'rpm' 'throttle_setting' 'current' 'power' 'motor_current' 'motor_power' 'propeller_efficiency' 'motor_efficiency' 'overall_efficiency'});
    end
    % Generate the full_throttle table for the aircraftdata struct
    tableSize = 41;
    T = zeros(tableSize,1); Q = zeros(tableSize,1); RPM = zeros(tableSize,1); P_tot = zeros(tableSize,1); P_motor = zeros(tableSize,1); eta_tot = zeros(tableSize,1);...
        eta_prop = zeros(tableSize,1); eta_motor = zeros(tableSize,1);I_motor = zeros(tableSize,1); U_motor = zeros(tableSize,1); I_ESC = zeros(tableSize,1);...
        U_ESC = zeros(tableSize,1); airspeed = zeros(tableSize,1);

    for i = 1:41
        airspeed(i) = i-1;
        [T(i), Q(i), RPM(i), P_tot(i), P_motor(i), ~ , eta_tot(i), eta_prop(i), eta_motor(i), ~, I_motor(i),...
            U_motor(i), I_ESC(i), U_ESC(i)] = interpolate_max_throttle_from_lookup(propulsion_lookup,ID_P,ID_M,ID_B,airspeed(i));
    end

    aircraftdata.Propulsion.full_throttle = table(airspeed,T,Q,RPM,U_motor./U_ESC, I_ESC, P_tot, I_motor, P_motor, eta_prop, eta_motor,eta_tot,'VariableNames',{'speed' 'thrust' 'torque' 'rpm' 'throttle_setting' 'current' 'power' 'motor_current' 'motor_power' 'propeller_efficiency' 'motor_efficiency' 'overall_efficiency'});
     % Fill the fields of the Propulsion part of the aircraftdata struct
    aircraftdata.Propulsion.battery = table2struct(Propulsion_Data_Struct.battery_table(ID_B,:));
    aircraftdata.Propulsion.motor = table2struct(Propulsion_Data_Struct.motor_table(ID_M,:));
    aircraftdata.Propulsion.propeller = Propulsion_Data_Struct.propeller_data(ID_P).data;
    aircraftdata.Propulsion.propeller.dynamic_table = Propulsion_Data_Struct.propeller_data(ID_P).dynamic;
    aircraftdata.Propulsion.propeller.static_table = Propulsion_Data_Struct.propeller_data(ID_P).static;

        % Error if the motor power exceeds the maximum rated power of the motor
    if aircraftdata.Propulsion.full_throttle.motor_power(1) > aircraftdata.Propulsion.motor.P_max
        error('Error202:maxPower','Maximum rated power of motor exceeded')
        % warning("Maximum rated power of motor exceeded")
    end
    
    % Error if ESC max current rating is exceeded
    if aircraftdata.Propulsion.full_throttle.motor_current(1) > 35
        error('Error203:excessESCcurrent','Motor current exceeds ESC continous current rating')
        % warning("Motor current exceeds ESC continous current rating")
    end

     % Error if battery discharge current exceeds the rated discharge current
    % of the battery
    if aircraftdata.Propulsion.full_throttle.current(1) > (aircraftdata.Propulsion.battery.C*aircraftdata.Propulsion.battery.Q)/1000
        error('Error204:excessBatteryDischargeCurrent','Battery current exceeds battery maximum discharge current')
        % warning("Battery current exceeds battery maximum discharge current")
    end
end