% Calculates the Liftcurve-Slope (C_l_alpha) from 2D-Polars
% Inputs: airspeed
% Outputs: Liftcurve-Slope for wing and tail airfoil
function [C_l_alpha_wing,C_l_alpha_tail] = liftCurveSlopeCalculation(aircraftdata,inputdata,airspeed)
%% Extract data from inputdata
density = inputdata.density;
viscosity = inputdata.viscosity;

%% Extract data from aircraftdata
mac_wing = aircraftdata.Geometry.wing.mean_aerodynamic_chord;
airfoilName_wing = aircraftdata.Geometry.wing.airfoil;

mac_tail = aircraftdata.Geometry.tail.mean_aerodynamic_chord;
airfoilName_tail = aircraftdata.Geometry.tail.airfoil;

chords = [mac_wing, mac_tail];
airfoilNames = {airfoilName_wing,airfoilName_tail};

C_l_alpha = zeros(1,2);

%% Calulate C_l_alpha for wing and tail airfoil
for k = 1:2
    airfoilName = airfoilNames{k};
    mac = chords(k);
    polars = inputdata.airfoil_data.(airfoilName).polar.Polars;
    Re_list = inputdata.airfoil_data.(airfoilName).polar.Re_list;

    C_l_alpha_list = zeros(1,length(Re_list));
    
    Re = density * airspeed * mac / viscosity;
    
    if Re > max(Re_list)
        error('Error102:ReNumber','Reynolds number too high')
    elseif Re < min(Re_list)
        error('Error102:ReNumber','Reynolds number too low')
    end
    
    for i = 1:length(Re_list)
        polar = polars{i};
        polar = polar(polar.CL > 0,:);
        for j = 1:size(polar,1)-1
            if (polar.CL(j+1) - polar.CL(j)) < 0
                break
            end
        end
        polar = polar(1:j,:);
        polar = polar(1:round(size(polar,1)/2),:);
        polar.alpha = polar.alpha * pi / 180;
        X = [ones(length(polar.alpha),1) polar.alpha];
        b = X\polar.CL;
        C_l_alpha_list(i) = b(2);
    end
    C_l_alpha(k) = interp1(Re_list,C_l_alpha_list,Re);
end
C_l_alpha_wing = C_l_alpha(1);
C_l_alpha_tail = C_l_alpha(2);
end