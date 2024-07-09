% Extracts maximum 2D lift Coeffcient for given Airfoil and speed
% Inputs: Reynold's Number (can be an array)
% Outputs: Maximum 2D-Lift Coefficient (same length as Re)
% TO DO:
%
function C_l_max = max2DLiftCoefficient(aircraftdata, inputdata, Re)

if any(strcmp(fieldnames(inputdata),'offsetdata'))
    offsetdata = inputdata.offsetdata;
    if any(strcmp(fieldnames(offsetdata),'Clmax'))
        offsetFactor = offsetdata.Clmax/100 + 1;
    else
        offsetFactor = 1;
    end
else
    offsetFactor = 1;
end

%% Extract data from input structs
airfoilName = aircraftdata.Geometry.wing.airfoil;
Re_list = inputdata.airfoil_data.(airfoilName).polar.Re_list;
polars = inputdata.airfoil_data.(airfoilName).polar.Polars;

%% Calculate max 2D-Lift Coefficient
C_l_max_list = zeros(1,length(Re_list));
for i = 1:length(Re_list)
    C_l_old = 0;
    for j = 1:length(polars{i}.CL)
        C_l = polars{i}.CL(j);
        if C_l_old > C_l && C_l > 0
            break;
        end
        C_l_old = C_l;
    end
    C_l_max_list(i) = C_l_old;
end

C_l_max = zeros(1,length(Re));
for j = 1:length(Re)
    C_l_max(j) = offsetFactor * interp1(Re_list,C_l_max_list,Re(j));
end
end