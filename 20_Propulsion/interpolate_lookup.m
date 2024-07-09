% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                            interpolate_lookup.m
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
% Interpolates the lookup table of a defined configuration to fit to a
% given speed
% -------------------------------------------------------------------------
% DATE          VERSION     AUTHOR              CHANGELOG
% 12/01/22    v1.0        Fabian Zach        initial release
% -------------------------------------------------------------------------
% Input variables:
% NAME          UNIT    DESCRIPTION
% propulsion_lookup     [-]         Struct of the lookup tables for the
%                                   propulsion system
% ID_P                  [-]         Propeller ID
% ID_M                  [-]         Motor ID
% ID_B                  [-]         Battery ID
% airspeed              [m s^-1]    The airspeed which the lookup table
%                                   should be interpolated for
% -------------------------------------------------------------------------
% Output variables:
% NAME                  UNIT        DESCRIPTION
% interp_lookup_table   [-]         Interpolated lookup table
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
function [interp_lookup_table] = interpolate_lookup(propulsion_lookup,ID_P,ID_M,ID_B,airspeed)
    % Do if the airspeed is not an integer value
    if airspeed > 0 && ~(mod(airspeed,1) == 0)
         % Get integer values of the airspeed next to the given airspeed
        airspeed1 = airspeed - mod(airspeed,1);
        airspeed2 = airspeed - mod(airspeed,1) + 1; 

        % Extract both airspeed tables from interpolation structures
        lookup1 = table2array(propulsion_lookup.prop(ID_P).motor(ID_M).battery(ID_B).airspeed(airspeed1 + 1).lookup_table);
        lookup2 = table2array(propulsion_lookup.prop(ID_P).motor(ID_M).battery(ID_B).airspeed(airspeed2 + 1).lookup_table);
        
        % Determine the sizes (length) of both extracted lookup tables
        lu1size = size(lookup1);
        lu2size = size(lookup2);

        % Set airspeed interpolation points
        airspeed_pts = [airspeed1 airspeed2];
        
        % The scope is required to handle the case if the adjacent lookup
        % tables are not of the same size
        scope = min(lu1size(1),lu2size(1));

        % Determine which lookup table is larger or if they are same size
        lu1_is_smaller = lu1size(1) < lu2size(1);
        same_size = lu1size(1) == lu2size(1);

        % Fill the first column of the new lookup table with RPM vlaues
        for i = 1:scope
            % Scope index is reached
            if i == scope

                % First lookup is smaller
                if lu1_is_smaller
                    % Interpolate the maxiumum RPM value
                    val = interp1(airspeed_pts,[lookup1(i,1) lookup2(i+1,1)],airspeed);
                    if val-i*100 > 0
                        % Fill the RPM in for the second to last row
                        % if it would not be filled
                        interp_mat(i,1) = i*100;
                        interp_mat(i+1,1) = val;
                    else
                        % Else just put in the RPM in the last row
                        interp_mat(i,1) = val;
                    end

                % Second lookup is smaller
                elseif ~lu1_is_smaller && ~same_size
                    % Interpolate the maxiumum RPM value
                    val = interp1(airspeed_pts,[lookup1(i+1,1) lookup2(i,1)],airspeed);
                    if val-i*100 > 0
                        % Fill the RPM in for the second to last row
                        % if it would not be filled
                        interp_mat(i,1) = i*100;
                        interp_mat(i+1,1) = val;
                    else
                        % Else just put in the RPM in the last row
                        interp_mat(i,1) = val;
                    end

                % Both lookup tables are of the same size
                elseif same_size                    
                    % Just put in the inerpolated RPM of the last row
                    interp_mat(i,1) = interp1(airspeed_pts,[lookup1(i,1) lookup2(i,1)],airspeed);
                end
            
            % If the scope is not reached just input the interpolated values of the two
            % lookup tables (this is unnecessary → change required)
            else
                interp_mat(i,1) = interp1(airspeed_pts,[lookup1(i,1) lookup2(i,1)],airspeed);
            end
        end
        
        % Fill the rest of the new lookup table with the interpolated
        % values
        for j = 2:17
            for i = 1:scope
                
                % End of the table not reached
                if i ~= scope
                    % Just interpolate the values between the two tables
                    % dependent on the airspeed
                    interp_mat(i,j) = interp1(airspeed_pts,[lookup1(i,j) lookup2(i,j)],airspeed);

                % The end is reached and the first lookup is smaller and scope is smaller than the length
                % of the bigger of the two tables
                elseif lu1_is_smaller && ~same_size && scope < max(lu1size(1),lu2size(1))
                    
                    % Interpolate the value for the last row of the new
                    % lookup table (from final value of LU1 and LU2)
                    interp_mat(end,j) = interp1(airspeed_pts,[lookup1(i,j) lookup2(end,j)],airspeed);
                    

                    % Interpolate the second to last value of the new
                    % lookup table
                    interp_mat(i,j) = interp1([interp_mat(i-2,1) interp_mat(i,1)], [interp_mat(end-2,j) interp_mat(end,j)], interp_mat(i,1));
                    interp_mat(i-1,1)
                 
                % The end is reached and the second lookup is smaller and scope is smaller than the length
                % of the bigger of the two tables
                elseif ~lu1_is_smaller && ~same_size && scope < max(lu1size(1),lu2size(1))

                    % Interpolate the value for the last row of the new
                    % lookup table (from final value of LU1 and LU2)
                    interp_mat(end,j) = interp1(airspeed_pts,[lookup1(end,j) lookup2(i,j)],airspeed);
                    
                    interp_mat(i,j) = interp1([interp_mat(i-2,1) interp_mat(i,1)], [interp_mat(end-2,j) interp_mat(end,j)], interp_mat(i,1));
                    interp_mat(i-1,1)
                 
                % Second lookup smaller and scope equals the length
                % of the smaller of the two tables
                elseif ~lu1_is_smaller && same_size && scope == min(lu1size(1),lu2size(1))
                    
                    interp_mat(end,j) = interp1(airspeed_pts,[lookup1(i,j) lookup2(i,j)],airspeed);
                end
            end
        end

        % Turn the matrix used for the interpolation into a table
        interp_lookup_table = array2table(interp_mat,"VariableNames",{'RPM' 'J' 'T' 'Q' 'Q_max' 'P_tot' 'eta_tot' 'eta_prop' 'P_shaft' 'P_motor' 'eta_motor' 'I_motor' 'V_motor' 'I_ESC' 'V_ESC' 'eta_ESC' 'eta_ESC_check'});

    % If the airspeed is an integer just extract the according lookup table
    % from the lookup struct
    elseif mod(airspeed,1) == 0
        interp_lookup_table = propulsion_lookup.prop(ID_P).motor(ID_M).battery(ID_B).airspeed(airspeed+1).lookup_table;
    end
    
end

