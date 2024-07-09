% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                             read_airfoil_polars.m
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
% Reads the saved Polars from XFOIL/XFLR and stores them in a struct
% -------------------------------------------------------------------------
% DATE          VERSION     AUTHOR              CHANGELOG
% 15/10/2019    v1.0        Thomas Seren        initial release
% -------------------------------------------------------------------------
% Input variables:
% NAME          UNIT    DESCRIPTION
% airfoil_name  -   name of used airfoil
% polar_path     -   path to the directory where the polar files are stored
% file_name_format - string definition of the polar file names BEHIND the airfoil name
% -------------------------------------------------------------------------
% Output variables:
% NAME          UNIT    DESCRIPTION
% Airfoil_Data_Struct  -  Struct with polar data
% -------------------------------------------------------------------------
% Limitations:  -
% Assumptions:  -
% -------------------------------------------------------------------------
% Example:      [Airfoil_Data_Struct, error_cell] = read_airfoil_polars('clarky', 'Polars', '_T1_Re%f_M0.00_N5.0_XtrTop25%_XtrBot25%.txt', cell(0,1))
% -------------------------------------------------------------------------
% Called functions:
% error_msg
% -------------------------------------------------------------------------
% Error Management
% ERROR MESSAGE                            DESCRIPTION
% -                                         -
% -------------------------------------------------------------------------
% KNOWN ISSUES:
% - 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [Airfoil_Data_Struct, error_cell] = read_airfoil_polars(airfoil_name, polar_path, file_name_format, error_cell)
    k = strfind(file_name_format,'%f');
    polar_listing = dir(fullfile(polar_path,[airfoil_name,'*',file_name_format(k+2:end)])); % find all files, that match the XFLR file name format
    if ~isempty(polar_listing)
        newest_polar_date = max([polar_listing.datenum]); % check newest filedate
        Airfoil_Data_file = fullfile(polar_path,[airfoil_name,'_Data.mat']); % .mat file name
        if exist(Airfoil_Data_file,'file')
            adf_listing = dir(Airfoil_Data_file);
            load(Airfoil_Data_file, 'Airfoil_Data_Struct')
            if newest_polar_date < adf_listing.datenum % check if .mat file is newer than newest polar
                return % return data from .mat file
            end
        else
            Airfoil_Data_Struct = struct('Airfoil',airfoil_name,'Re_list',[],'Polars',[]);

        end
        Re_list = Airfoil_Data_Struct.Re_list; % list of all existing Reynolds numbers
        Polars =  Airfoil_Data_Struct.Polars; % all polars
        if isempty(Polars)
            Polars = cell(0,1);
        end
        for i_file = 1:numel(polar_listing) % go through all polar data files
            curname = polar_listing(i_file).name;
            cur_Re = sscanf(curname,[airfoil_name,file_name_format]) * 1e6;
            if ~ismember(cur_Re, Re_list) % if polar does not exist already load data and put it in the polar table cell
                curfile = fullfile(polar_listing(i_file).folder, polar_listing(i_file).name);
                Re_list(end+1) = cur_Re;
                Input = dlmread(curfile,'%[ ]',11,0);
                Input((Input(:,3)==0),:)=[];
                [~, CL_max] = max(Input(:,2)); %find the position of CL_max
                [~, CL_min] = min(Input(:,2));% find the position of CL_min
                Input = Input(CL_min:CL_max,:); %reduce to the range between CL_min and CL_max
                Input = make_monotonic(Input,2); % sort out values that are not strictily monotonic increasing
                Polar_Table = array2table(Input, 'VariableNames', {'alpha','CL','CD','CDp','Cm','TopXtr','BotXtr','Cpmin','Chinge','XCp'});
                Polars{end+1} = Polar_Table;
            end
        end
        [Re_list, I] = sort(Re_list); % sort the reynlods number list
        Polars = Polars(I); % so
        Airfoil_Data_Struct.Polars = Polars;
        Airfoil_Data_Struct.Re_list = Re_list;
        save(Airfoil_Data_file, 'Airfoil_Data_Struct')
    else
        msg = ['No polar data for airfoil ', airfoil_name,' found!'];
        error_cell = error_msg(msg, error_cell);
        Airfoil_Data_Struct = struct();
    end
end
function matrix = make_monotonic(matrix, column)
    if ~exist('column','var')
        column = 1;
    end
    sel = monotonic_selection(matrix(:,column));
    while ~all(~sel)
        matrix(sel,:) = [];
        sel = monotonic_selection(matrix(:,column));
    end
end
function sel = monotonic_selection(vector)
    vector = reshape(vector,1,[]);
    sel = logical([0, ~(vector(1:end-1) < vector(2:end))]);
end