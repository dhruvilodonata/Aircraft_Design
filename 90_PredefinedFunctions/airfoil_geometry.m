% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                             airfoil_geometry.m
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
% Calculates the relative crosssection, circumference of the airfoil
% Calculates the height of the wing spar
% -------------------------------------------------------------------------
% DATE          VERSION     AUTHOR              CHANGELOG
% 15/10/2019    v1.0        Thomas Seren        initial release
% -------------------------------------------------------------------------
% Input variables:
% NAME          UNIT    DESCRIPTION
% airfoil_filename  -   file name of used airfoil
% spar_Position     -   relative position of the spar
% -------------------------------------------------------------------------
% Output variables:
% NAME          UNIT    DESCRIPTION

% -------------------------------------------------------------------------
% Limitations:  -
% Assumptions:  -
% -------------------------------------------------------------------------
% Example:      [geo_struct, error_cell] = airfoil_geometry('clarky', 0.25, cell(0,1))
% -------------------------------------------------------------------------
% Called functions:
% InterX.m from matlab fileexchange https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections
% error_msg.m
% -------------------------------------------------------------------------
% Error Management
% ERROR MESSAGE                            DESCRIPTION
% 'No airfoil geometry for airfoil found!']     there is no airfoil.dat at the specified location
% -------------------------------------------------------------------------
% KNOWN ISSUES:
% - airfoil coordinates must begin at the upper end of the airfoil, run until the leading edge is reached and then back to the lower end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [geo_struct, error_cell] = airfoil_geometry(airfoil_name, airfoil_directory, spar_position, error_cell)
    airfoil_filename = fullfile(airfoil_directory, [airfoil_name,'.dat']);
    if exist(airfoil_filename,'file')
        fid_airfoil = fopen(airfoil_filename); % open textfile
        input_array_airfoil_data = textscan(fid_airfoil, '%f %f', 1000, 'headerlines', 1); %read text
        airfoil_coordinate_matrix = cell2mat(input_array_airfoil_data); %create array from input
        
        spar_coordinates = InterX(airfoil_coordinate_matrix.',[spar_position * [1,1]; -1 1]); % calculate the spar top and bottom coordinats
        
        % relative airfoil thickness spar to airfoil length
        relative_spar_thickness = spar_coordinates(2,2) - spar_coordinates(2,1);
        
        % relative airfoil crosssection to airfoil length^2
        relative_crosssection = (airfoil_coordinate_matrix(1:end-1,1)-airfoil_coordinate_matrix(2:end,1)).'*(airfoil_coordinate_matrix(1:end-1,2)+airfoil_coordinate_matrix(2:end,2))/2;
        
        % relative x position of the CG of the crosssection
        crosssection_CG = 1/6/relative_crosssection * (airfoil_coordinate_matrix(1:end-1,1)+airfoil_coordinate_matrix(2:end,1)).'*(airfoil_coordinate_matrix(1:end-1,1).*airfoil_coordinate_matrix(2:end,2)-airfoil_coordinate_matrix(2:end,1).*airfoil_coordinate_matrix(1:end-1,2));
        
        relative_circumference_vector = (((airfoil_coordinate_matrix(1:end-1,1) - airfoil_coordinate_matrix(2:end,1)).^2+(airfoil_coordinate_matrix(1:end-1,2) - airfoil_coordinate_matrix(2:end,2)).^2)).^.5;
        
        % relative airfoil circumference to airfoil length
        relative_circumference = sum(relative_circumference_vector); 

        xcg_pos = (airfoil_coordinate_matrix(1:end-1,1) + airfoil_coordinate_matrix(2:end,1))/2;
        
        % relative x position of the CG of the skin
        skin_CG = (xcg_pos.'*relative_circumference_vector)/sum(relative_circumference);
        
        
        geo_struct = struct(...
            'relative_circumference' ,relative_circumference, ...
            'relative_crosssection' ,relative_crosssection, ...
            'relative_spar_thickness' ,relative_spar_thickness, ...
            'skin_CG' ,skin_CG, ...
            'crosssection_CG' ,crosssection_CG ...
            );
    else
        geo_struct = struct();
        msg = ['No airfoil geometry for airfoil ', airfoil_name, ' found!'];
        error_cell = error_msg(msg, error_cell);
    end
end

