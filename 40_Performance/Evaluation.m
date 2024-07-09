clc
clear
% Function to evaluate and illustrate the results for aircraft configurations
currentFolder = pwd;
contents = dir(fullfile(currentFolder,'70_Results'));
filename = '202201232221_combination_table.mat';
loadedFile = load(fullfile(currentFolder,'70_Results',filename));
combination_table = loadedFile.combination_table;

fieldsToPlot = {'wing_airfoil', 'wing_S_ref', 'wing_aspectRatio', 'propeller_ID', 'motor_ID', 'battery_ID'};

motor_id_array = unique(combination_table.motor_ID)';
propeller_id_array = unique(combination_table.propeller_ID)';
motor_propeller_array = combvec(motor_id_array,propeller_id_array)';

for i = 1:size(motor_propeller_array,1)
    motor_i = combination_table(combination_table.motor_ID == motor_propeller_array(i,1),:);
    combination_i = motor_i(motor_i.propeller_ID == motor_propeller_array(i,2),:);
    errors = combination_i.errorID;
end

for i = 1:length(fieldsToPlot)
    %subplot(subplotRows,subplotColoumns,i)
    figure
    if strcmp('wing_airfoil',fieldsToPlot{i})
        airfoil = combination_table.wing_airfoil;
        airfoil_list = unique(airfoil);
        xdata = zeros(length(airfoil),1);
        for j = 1:length(airfoil_list)
            xdata = xdata + j * strcmp(airfoil,airfoil_list(j));
        end
        %xdata = airfoil();
    else
        xdata = combination_table.(fieldsToPlot{i});
    end
    plot(xdata,combination_table.points,'*')
    xlabel(strrep(fieldsToPlot{i},'_',' '))
    ylabel('Points')
end