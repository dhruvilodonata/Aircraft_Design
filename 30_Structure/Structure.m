function [ aircraftdata, error_cell ] = Structure( aircraftdata, inputdata, error_cell)
%% Set constant values
sigma_m = 1500e6; % [Pa]
G = 670e6; %shear modulus [Pa]
airfoil_skin_thickness = 0.8e-3; % [m]
tau_allowable = 15.5e6; % [Pa]

%% Extract data from aircraft data struct
H_airfoil= aircraftdata.Geometry.wing.thickness;
b = aircraftdata.Geometry.wing.span;
taper_ratio = aircraftdata.Geometry.wing.taper_ratio;
mac_wing = aircraftdata.Geometry.wing.mean_aerodynamic_chord;
airfoil_wing = aircraftdata.Geometry.wing.airfoil;
y = aircraftdata.Aerodynamics.wing_loading.v_max_CL_max.load_table.y_position;
wing_load_data = aircraftdata.Aerodynamics.wing_loading;

%% Extract data from inputdata
wing_relative_crosssection = inputdata.airfoil_data.(airfoil_wing).geo.relative_crosssection;
circ_relative = inputdata.airfoil_data.(airfoil_wing).geo.relative_circumference; 

chord_root = 2 * mac_wing / (1 + taper_ratio);
chord_array = chord_root - ( chord_root - chord_root* taper_ratio)/(b/2)*y;

%% Bending Moment across the wing
cases = {'stall', 'v_max_CL_max'};
for index = 1:2
    loadCase = cases{index};
    C_l = wing_load_data.(loadCase).load_table.C_l;
    q = wing_load_data.(loadCase).dynamic_pressure;
    if isempty(C_l)
        continue
    end
    shear_force = zeros(length(y),1);
    bending_moment = zeros(length(y),1);
    
    lift_line_force = q*C_l.*chord_array;
    a = [0;y];
    for i=1:length(y)
        for j = 1:i
        shear_force(i) = shear_force(i) + lift_line_force(j)*(a(j+1)-a(j));
        end
    end
    shear_force = shear_force - shear_force(length(shear_force));
    for i=1:length(bending_moment)
        for j = 1:i
        bending_moment(i) = bending_moment(i) + shear_force(j)*(a(j+1)-a(j));
        end
    end
    bending_moment = bending_moment - bending_moment(length(shear_force));
    
    lift_line_distribution.(loadCase) = lift_line_force;
    shear_distribution.(loadCase) = shear_force;
    bending_moment_distribution.(loadCase) = bending_moment;
end
%% Bending Stress
Z_distance = H_airfoil/2; %max stress distance from neutral point
bending_moment_max = max(bending_moment_distribution.v_max_CL_max);
radius = [0.0005 0.00075 0.001];
J_inertia = 2 * ((pi*radius.^2).*(H_airfoil/2 - radius).^2);
sigma = (bending_moment_max*Z_distance./J_inertia);% Pa
Safety_factor = sigma_m./sigma;

if any(Safety_factor > 1.5)
    index = find(Safety_factor > 1.5);
    index = index(1);
    spar_diameter = radius(index) * 2;
else
    error('Error303:safetyFactor','Critical Safety Factor of 1.5 could not be reached for bending moment')
end

%% Shear Force
Q_max = max(abs(shear_distribution.v_max_CL_max));
thickness_sparweb = (1.5*Q_max)/(tau_allowable*(H_airfoil-2*spar_diameter));
%A_sparweb = thickness_sparweb*(H_airfoil-2*spar_diameter);%m^2

%% Wing Torsion
wing_crosssection = wing_relative_crosssection * chord_array.^2;
circumference = circ_relative * chord_array;

% maxspeed case
q = wing_load_data.v_max_CL_max.dynamic_pressure;
C_m = wing_load_data.v_max_CL_max.load_table.C_m;
torsional_force = q*chord_array.^2.*C_m; %Torsion N/m^2
torsional_moment_distribution = integral(torsional_force,y);
I_T = 4*wing_crosssection.^2 ./ circumference * airfoil_skin_thickness;
phi_dot = torsional_moment_distribution ./ (I_T * G);
twist = flip(integral(flip(phi_dot),flip(y)));
twist_max_deg = max(twist) * 180 / pi;
twist_SF = 4 / abs(twist_max_deg);

if twist_SF < 1.5
    error('Error304:twistAngle','Twist Angle Exceeds Critical Angle of 4 degrees')
end

%% Write Results to Aircraft Data Struct
aircraftdata.Structure.maximum_wing_load_table.y_position = y;
aircraftdata.Structure.maximum_wing_load_table.q = lift_line_distribution.v_max_CL_max;
aircraftdata.Structure.maximum_wing_load_table.Q = shear_distribution.v_max_CL_max;
aircraftdata.Structure.maximum_wing_load_table.M_b = bending_moment_distribution.v_max_CL_max;
aircraftdata.Structure.maximum_wing_load_table.M_t = torsional_moment_distribution;

aircraftdata.Structure.spar_diameter = spar_diameter;     % [m]
aircraftdata.Structure.spar_web_thickness = thickness_sparweb; % [m]
aircraftdata.Structure.wing_twist = twist_max_deg;       % [deg]
end