% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                             read_propulsion_data.m
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
% Reads the propulsion element data from excel and saves them in a struct
% -------------------------------------------------------------------------
% DATE          VERSION     AUTHOR              CHANGELOG
% 15/10/2019    v1.0        Thomas Seren        initial release
% -------------------------------------------------------------------------
% Input variables:
% NAME          UNIT    DESCRIPTION
% directory  -   directory, where the propulsion data is stored
% error_cell     -   cell to collect all error messages
% -------------------------------------------------------------------------
% Output variables:
% NAME          UNIT    DESCRIPTION
% Propulsion_Data_Struct  -  Struct with propulsion data
% error_cell     -   cell to collect all error messages
% -------------------------------------------------------------------------
% Limitations:  -
% Assumptions:  -
% -------------------------------------------------------------------------
% Example:      [Propulsion_Data_Struct, error_cell] = read_propulsion_data('../80_Databses/Propulsion', cell(0,1))
% -------------------------------------------------------------------------
% Called functions:
% -
% -------------------------------------------------------------------------
% Error Management
% ERROR MESSAGE                            DESCRIPTION
% 'File not found!'     the specified files can't be found - check if the
%                       directory is correct
% -------------------------------------------------------------------------
% KNOWN ISSUES:
% - 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%
function [Propulsion_Data_Struct, error_cell] = read_propulsion_data(directory, error_cell)
    %%
    [Propulsion_Data_Struct.motor_table, error_cell] = overviewinput('Motor_List', directory, error_cell);
    [Propulsion_Data_Struct.battery_table, error_cell] = overviewinput('Battery_List', directory, error_cell);
    [Propulsion_Data_Struct.propeller_table, error_cell] = overviewinput('Propeller_List', directory, error_cell);
    [Propulsion_Data_Struct.propeller_data, error_cell] = propellerdata_input(Propulsion_Data_Struct.propeller_table, directory, error_cell);

end
function [data_table, error_cell] = overviewinput(filename, directory, error_cell)
    %%
    excelpath = fullfile(directory, [filename,'.xlsx']); % path to excel table
    if ~exist(excelpath, 'file') % check if excel file can be found
        msg = ['File ', excelpath,' not found!'];
        error_cell = error_msg(error_cell,msg);
        data_table = table();
        return
    end
    matpath = fullfile(directory, [filename,'.mat']); % path to matlab workspace (.mat) file
    if exist(matpath, 'file') % check if .mat file exists
        xlsdir = dir(excelpath); % get file information of excel file
        matdir = dir(matpath); % get file information of .mat file
        if xlsdir.datenum < matdir.datenum % check if there had been any changes in excel file (i.e. it is newer)
            load(matpath, 'data_table'); % if .mat file is newer than excel file load .mat file --> is much faster
            return
        end
    end
    [~, ~, raw] = xlsread(excelpath,'overview'); %load excel file
    data_table = cell2table(raw(4:end,2:end),'VariableNames',raw(3,2:end)); % create data table from excel file
    save(matpath, 'data_table'); % save updated/created .mat file
    
end
function [propeller_data, error_cell] = propellerdata_input(propeller_table, directory, error_cell)
%%
    excelpath = fullfile(directory,'Propeller_List.xlsx'); % path to excel table
    if ~exist(excelpath, 'file') % check if excel file can be found
        msg = ['File ', excelpath,' not found!'];
        error_cell = error_msg(error_cell,msg);
        propeller_data = struct();
        return
    end
    matpath = fullfile(directory, 'Propeller_Polars.mat'); % path to matlab workspace (.mat) file
    if exist(matpath, 'file') % check if .mat file exists
        xlsdir = dir(excelpath); % get file information of excel file
        matdir = dir(matpath); % get file information of .mat file
        if xlsdir.datenum < matdir.datenum % check if there had been any changes in excel file (i.e. it is newer)
            load(matpath, 'propeller_data');  % if .mat file is newer than excel file load .mat file --> is much faster
            return
        end
    end
    n_propeller = size(propeller_table,1); % number of proepller
    propeller_data = repmat(struct('data',table,'static',table,'dynamic',table),n_propeller,1); % create propeller data struct
    for i_propeller = 1:n_propeller
        [~, ~, raw] = xlsread(excelpath,num2str(i_propeller,'ID%i')); % read propeller polar data
        propeller_data(i_propeller).dynamic = cell2table(raw(2:end,1:4),'VariableNames',raw(1,1:4)); % create dynamic table
        propeller_data(i_propeller).static = cell2table(raw(2:end,end-2:end),'VariableNames',raw(1,end-2:end)); % create static table
        propeller_data(i_propeller).static(isnan(table2array(propeller_data(i_propeller).static(:,1))),:) = []; %clean static table
        propeller_data(i_propeller).data = table2struct(propeller_table(i_propeller,:)); % add propeller information
    end
    save(matpath, 'propeller_data'); % save updated/created .mat file

end