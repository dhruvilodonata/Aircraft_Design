% Calculates lift distribution, induced drag coeffcient and angle of attack from given C_L in AVL
% TO DO:
function [y,C_l,C_m,l_C_l,C_D_i,AoA,stall,position_NP] = AVL_Call(aircraftdata,inputdata,airspeed,C_L)
%% Extract information from input structs
viscosity = inputdata.viscosity;
density = inputdata.density;
output = inputdata.AVLOutput;
identifier = aircraftdata.id;

mac_wing = aircraftdata.Geometry.wing.mean_aerodynamic_chord;
taperRatio_wing = aircraftdata.Geometry.wing.taper_ratio;
b_wing = aircraftdata.Geometry.wing.span;

if output
    outputString = '';
else
    outputString = ' > NUL';
end
currentPath = pwd;
if ispc
    avlFolderPath = fullfile(currentPath,'10_Aerodynamics\AVL');
    avlPath = fullfile(avlFolderPath,'avl.exe');
    systemCommand = ['"',avlPath,'" < "',fullfile(avlFolderPath,['AVLScript',num2str(identifier),'.txt']),'"',outputString];
else
    avlFolderPath = fullfile(currentPath,'10_Aerodynamics/AVL');
    avlPath = fullfile(avlFolderPath,'avl3.35');
    systemCommand = [strrep(avlPath,' ','\ '),' < ', strrep(fullfile(avlFolderPath,['AVLScript',num2str(identifier),'.txt']),' ','\ '),outputString];
end
%% Create AVL command file
index = 0;
while true
    index = index + 1;
    executefile = fopen(fullfile(avlFolderPath,['AVLScript',num2str(identifier),'.txt']),'w');
    if ~(executefile < 0)
        break
    elseif index > 5
        error('Error999:FileID',['Could not open AVLScript file for combination ',identifier])
    end
end
fileinput = {...
    'plop',...
    'g',... % Deactivating graphical output
    '',...
    ['LOAD ',fullfile(avlFolderPath, ['ourOwnPlane',num2str(identifier),'.avl'])],...
    'oper',...
    'm',...
    'd',...
    num2str(density),...
    'v',...
    num2str(airspeed),...
    '',...
    'd1 pm 0',...
    ['a c ',num2str(C_L)],...
    'x',...
    'FT',...
    ['totalForces',num2str(identifier),'.ft'],...
    'FS',...
    ['stripForces',num2str(identifier),'.fs'],...
    'ST',...
    ['stabilityDerivatives',num2str(identifier),'.st']...
};
if ispc
    fprintf(executefile,strrep(strjoin(fileinput,'\n'),'\','\\'));
else
    fprintf(executefile,strjoin(fileinput,'\n'));
end
index = 1;
while true
    status = fclose(executefile);
    if status == 0
        break
    elseif index > 5
        error('Error999:FileID',['Could not close AVLScript file for combination ',identifier])
    end
    index = index + 1;
end
%% Run AVL
wingIncidenceAngle = 0;
i = 0;
while true
    i = i + 1;
    createAVLFile(aircraftdata,inputdata,wingIncidenceAngle,identifier,airspeed);
    if output
        system(systemCommand);
    else
        [~,~] = system(systemCommand);
    end
    C = readcell(['totalForces',num2str(identifier),'.ft'],'FileType','text','Delimiter',' ',...
    'ConsecutiveDelimitersRule','join','LeadingDelimitersRule','ignore','Range','30:30');
    deflection = C{1,3}; %read AVL file
    wingIncidenceAngle = wingIncidenceAngle + deflection;
    if abs(deflection) > 0.5
        delete(['stripForces',num2str(identifier),'.fs'],['totalForces',num2str(identifier),'.ft'],['stabilityDerivatives',num2str(identifier),'.st']);
    elseif i > 10
        error('Error101:TrimError','Did not converge. A/C not trimmable');
    else
        break
    end
end
%% Read Lift Distribution from AVL
C = readcell(['stripForces',num2str(identifier),'.fs'],'FileType','text','Delimiter',' ',...
    'ConsecutiveDelimitersRule','join','LeadingDelimitersRule','ignore','Range','21:40');
y = cell2mat(C(:,2));
C_l = cell2mat(C(:,7));
C_m = cell2mat(C(:,11));
l_C_l = zeros(length(C_l),1);%cell2mat(C(:,13));
%% Read Induced Drag from AVL
C = readcell(['totalForces',num2str(identifier),'.ft'],'FileType','text','Delimiter',' ',...
    'ConsecutiveDelimitersRule','join','LeadingDelimitersRule','ignore','Range','26:26');
C_D_i = C{1,6};
%% Read Angle of attack from AVL
C = readcell(['totalForces',num2str(identifier),'.ft'],'FileType','text','Delimiter',' ',...
    'ConsecutiveDelimitersRule','join','LeadingDelimitersRule','ignore','Range','16:16');
AoA = cell2mat(C(1,3));
%% Read Neutral Point Position from AVL
C = readcell(['stabilityDerivatives',num2str(identifier),'.st'],'FileType','text','Delimiter',' ',...
    'ConsecutiveDelimitersRule','join','LeadingDelimitersRule','ignore','Range','64:64');
position_NP = round(cell2mat(C(1,5)),2);
%% Delete created Files
delete(fullfile(avlFolderPath,['ourOwnPlane',num2str(identifier),'.avl']),fullfile(avlFolderPath,['AVLScript',num2str(identifier),'.txt']));
delete(['stripForces',num2str(identifier),'.fs'],['totalForces',num2str(identifier),'.ft'],['stabilityDerivatives',num2str(identifier),'.st']);
%% Check if stall occures
chord_root_wing = 2 * mac_wing / (1 + taperRatio_wing);
chord = chord_root_wing - (chord_root_wing - chord_root_wing * taperRatio_wing)/(b_wing/2)*y;
Re = density * airspeed * chord / viscosity;
C_l_max = max2DLiftCoefficient(aircraftdata, inputdata, Re);
stall = any(C_l' > C_l_max);
end