function [aircraftdata, error_cell] = Precalculation(aircraftdata,inputdata,error_cell)
if any(strcmp(fieldnames(inputdata),'offsetdata'))
    offsetdata = inputdata.offsetdata;
    if any(strcmp(fieldnames(offsetdata),'S_ref'))
    offsetFactor_S_ref = offsetdata.S_ref/100 + 1;
    else
        offsetFactor_S_ref = 1;
    end
    if any(strcmp(fieldnames(offsetdata),'AR'))
        offsetFactor_AR = offsetdata.AR/100 + 1;
    else
        offsetFactor_AR = 1;
    end
else
    offsetFactor_S_ref = 1;
    offsetFactor_AR = 1;
end

%% Extract data from aircraft data struct
S_ref = offsetFactor_S_ref * aircraftdata.Configuration.design_parameter.wing_S_ref;
AR_wing = offsetFactor_AR * aircraftdata.Configuration.design_parameter.wing_aspectRatio;
taperRatio_wing = aircraftdata.Configuration.fixed_parameter.wing_taperRatio;

AR_tail = aircraftdata.Configuration.fixed_parameter.tail_aspectRatio;
taperRatio_tail = aircraftdata.Configuration.fixed_parameter.tail_taperRatio;

airfoil_wing = aircraftdata.Configuration.design_parameter.wing_airfoil;
airfoil_wing_relative_crosssection = inputdata.airfoil_data.(airfoil_wing).geo.relative_crosssection;
airfoil_wing_relative_circumference = inputdata.airfoil_data.(airfoil_wing).geo.relative_circumference;
airfoil_wing_relative_spar_thickness = inputdata.airfoil_data.(airfoil_wing).geo.relative_spar_thickness;

airfoil_tail = aircraftdata.Configuration.fixed_parameter.tail_airfoil;
dihedral = aircraftdata.Configuration.fixed_parameter.dihedral;
dihedral_position = aircraftdata.Configuration.fixed_parameter.dihedral_position;

%% Calculations
b_wing = sqrt(AR_wing * S_ref);
mac_wing = S_ref / b_wing;
chord_root_wing = 2 * mac_wing / (1 + taperRatio_wing);

wing_seperation_y = [0,b_wing/2];
wing_seperation_chord = chord_root_wing - (chord_root_wing - chord_root_wing * taperRatio_wing)/(b_wing/2)*wing_seperation_y;
wing_seperation_crosssections = airfoil_wing_relative_crosssection * wing_seperation_chord.^2;     % per wing
wing_seperation_crosssections_sum = sum(wing_seperation_crosssections); % per wing

volume_wing = airfoil_wing_relative_crosssection * (0.5 * b_wing * chord_root_wing^2 - chord_root_wing^2 * (1 - taperRatio_wing) * b_wing/2 ...
                + chord_root_wing^2 * (1 - taperRatio_wing)^2 * b_wing / 6); % Is this correct?
S_wet_wing = airfoil_wing_relative_circumference * mac_wing * b_wing + 2 * wing_seperation_crosssections_sum;

wing_thickness = airfoil_wing_relative_spar_thickness * chord_root_wing;

%% Write results to aircraft data struct
aircraftdata.Geometry.wing.aspect_ratio = AR_wing;
aircraftdata.Geometry.wing.projected_surface = S_ref;
aircraftdata.Geometry.wing.span = b_wing;
aircraftdata.Geometry.wing.mean_aerodynamic_chord = mac_wing;
aircraftdata.Geometry.wing.root_aerodynamic_chord = chord_root_wing;
aircraftdata.Geometry.wing.taper_ratio = taperRatio_wing;
aircraftdata.Geometry.wing.volume = volume_wing;
aircraftdata.Geometry.wing.wetted_surface = S_wet_wing;
aircraftdata.Geometry.wing.thickness = wing_thickness;
aircraftdata.Geometry.wing.airfoil = airfoil_wing;
aircraftdata.Geometry.wing.dihedral = dihedral;
aircraftdata.Geometry.wing.dihedral_position = dihedral_position;

aircraftdata.Geometry.tail.airfoil = airfoil_tail;
aircraftdata.Geometry.tail.aspect_ratio = AR_tail;
aircraftdata.Geometry.tail.taper_ratio = taperRatio_tail;
end

