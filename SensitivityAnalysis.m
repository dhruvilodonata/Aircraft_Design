%% Preparation
clear
clc
load(fullfile('70_Results','Final_Configuration.mat'))
inputdata_default = inputdata;
points_benchmark = aircraftdata.Performance.points;
error_cell = cell(1,1);
aircraftdata.Performance.stall_speed = NaN;
sensitivityParameters = {'mass','drag','Clmax','S_ref','AR','power'};
offsets = -10:1:10;
T = array2table(zeros(length(offsets),length(sensitivityParameters)),'VariableNames',sensitivityParameters);
T.offset = offsets';
T = movevars(T,'offset','Before','mass');
%% Calculation
for i = 1:length(sensitivityParameters)
    for j = 1:length(offsets)
        if offsets(j) ~= 0
            inputdata = inputdata_default;
            inputdata.offsetdata.(sensitivityParameters{i}) = offsets(j);
            [aircraftdata_new,error_cell] = objective_function(aircraftdata,inputdata,error_cell);
            points = aircraftdata_new.Performance.points;
            T.(sensitivityParameters{i})(j) = (points / points_benchmark - 1) * 100;
            if aircraftdata_new.Performance.stall_speed > 9
                T.(sensitivityParameters{i})(j) = T.(sensitivityParameters{i})(j)+1000;
            end
        else
            T.(sensitivityParameters{i})(j) = 0;
        end
    end
end
%% Adjust Table for stallSpeed
for i  = 1:length(sensitivityParameters)
    stall = T.(sensitivityParameters{i}) > 900;
    if any(stall)
        T.([sensitivityParameters{i},'_stall']) = T.(sensitivityParameters{i})-1000;
        T.([sensitivityParameters{i},'_stall'])(~stall) = NaN;
        T.(sensitivityParameters{i})(stall) = NaN;
    end
end
%% Plots
colorCell = {'#0072BD','#D95319','#A2142F','#7E2F8E','#77AC30','#4DBEEE'};
legendEntries = {'Mass','Drag','C_{l,max}','S_{ref}','Aspect Ratio','Power'};
figure
hold on
h = zeros(length(sensitivityParameters),1);
for i = 1:length(sensitivityParameters)
    h(i) = plot(T.offset,T.(sensitivityParameters{i}),'Color',colorCell{i},'DisplayName',legendEntries{i});
    if any(strcmp(T.Properties.VariableNames,[sensitivityParameters{i},'_stall']))
        plot(T.offset,T.([sensitivityParameters{i},'_stall']),':','Color',colorCell{i},'DisplayName',legendEntries{i});
    end
end
legend(h)
xlabel('Deviation of Parameter (in %)')
ylabel('Deviation in Points achieved (in %)')
grid()