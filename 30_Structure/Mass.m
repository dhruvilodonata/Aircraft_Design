function [ aircraftdata, error_cell ] = Mass( aircraftdata, inputdata, error_cell)

if any(strcmp(fieldnames(inputdata),'offsetdata'))
    offsetdata = inputdata.offsetdata;
    if any(strcmp(fieldnames(offsetdata),'mass'))
        offsetFactor = offsetdata.mass/100 + 1;
    else
        offsetFactor = 1;
    end
else
    offsetFactor = 1;
end

%% Extract data
length_RC = 175e-3;
position_tail = aircraftdata.Geometry.tail.x_position;

S_wet_wing = aircraftdata.Geometry.wing.wetted_surface;
V_wing = aircraftdata.Geometry.wing.volume;
S_wet_TP = aircraftdata.Geometry.tail.wetted_surface;
V_TP = aircraftdata.Geometry.tail.volume;

b_wing = aircraftdata.Geometry.wing.span;
taper_wing = aircraftdata.Geometry.wing.taper_ratio;
spar_diameter = aircraftdata.Structure.spar_diameter;
mac_wing = aircraftdata.Geometry.wing.mean_aerodynamic_chord;
wing_thickness = aircraftdata.Geometry.wing.thickness;
taper_ratio_wing = aircraftdata.Geometry.wing.taper_ratio;
S_ref = aircraftdata.Geometry.wing.projected_surface;
airfoilName = aircraftdata.Geometry.wing.airfoil;
thickness_airfoil_relative = inputdata.airfoil_data.(airfoilName).geo.relative_spar_thickness;
thickness_web = aircraftdata.Structure.spar_web_thickness;
skin_CG = inputdata.airfoil_data.(airfoilName).geo.skin_CG;

position_cg_target = aircraftdata.Aerodynamics.target_position_CG;

mass_propeller = 1e-3 * aircraftdata.Propulsion.propeller.m;
mass_motor = 1e-3 * aircraftdata.Propulsion.motor.m;
mass_battery = 1e-3 * aircraftdata.Propulsion.battery.m;

% Part Names
parts = {'ESC wire';'Fuselage';'Tailplane Connector';'Wing Separation';'Main Landing Gear';...
    'Wing Pylon'; 'Wing'; 'Spar Caps'; 'Tailplane';'Propeller';'Motor';'Battery';...
    'Tailplane Adaptor'; 'Main Landing Gear Adaptor';'Wheel';'Motor mount';'RC-Box';...
    'Wing Actuators';'Wing Connector';'Spar web'};

position_fuselage_tip = -0.2;
position_RC_offset = 0;
i = 1;
j = 1;
while true
    % Part Position
    integral_circumference_over_wing =...
    S_ref / (1+taper_wing)*(1-(1-taper_wing)/2);

    position_wing = (skin_CG-0.25) * (2*S_ref / (1+taper_wing)/b_wing)^2*...
        (b_wing/2-(1-taper_wing)*b_wing/2+(1-taper_wing)^2*b_wing/6)...
        /integral_circumference_over_wing;

    positions = zeros(length(parts),1);
    positions(strcmp(parts,'Wing'),1) = position_wing;
    positions(strcmp(parts,'Tailplane'),1) = position_tail;
    positions(strcmp(parts,'Tailplane Connector'),1) = position_tail;
    positions(strcmp(parts,'Tailplane Adaptor'),1) = position_tail;
    positions(strcmp(parts,'Wing Actuators'),1) = mac_wing * 0.25;
    positions(strcmp(parts,'Wing Pylon'),1) = mac_wing * 0.25;
    positions(strcmp(parts,'Fuselage'),1) = 0.5 * (position_fuselage_tip + position_tail);
    positions(strcmp(parts,'Propeller'),1) = position_fuselage_tip - 75e-3;
    positions(strcmp(parts,'Motor'),1) = position_fuselage_tip - 52.5e-3;
    positions(strcmp(parts,'Motor mount'),1) = position_fuselage_tip;
    positions(strcmp(parts,'Battery'),1) = position_fuselage_tip + length_RC/2;
    positions(strcmp(parts,'RC-Box'),1) = position_fuselage_tip + 20e-3 + length_RC/2 + position_RC_offset;
    positions(strcmp(parts,'ESC wire'),1) = 0.5*(positions(strcmp(parts,'Motor'),1) + positions(strcmp(parts,'RC-Box'),1));
    
    % Part Quantities
    quantities = ones(length(parts),1);
    quantities(strcmp(parts,'Tailplane Connector'),1) = 2;
    quantities(strcmp(parts,'Wing Separation'),1) = 2 * (ceil(b_wing / 2 / 0.58)-1);
    quantities(strcmp(parts,'ESC Wire'),1) = 0.25;
        %abs(positions(strcmp(parts,'RC-Box'),1) - positions(strcmp(parts,'Motor'),1));
    quantities(strcmp(parts,'Fuselage'),1) = 1.122;...
        %abs(positions(strcmp(parts,'Tailplane'),1) - positions(strcmp(parts,'Motor'),1));
    quantities(strcmp(parts,'Tailplane Adaptor'),1) = 1;
    quantities(strcmp(parts,'Wheel'),1) = 2;
    quantities(strcmp(parts,'Spar Caps'),1) = b_wing * 2;
    
    % Part Mass
    if isnan(thickness_web)
        thickness_web = 0;
    end
    chord_root = 2 * mac_wing / (1 + taper_ratio_wing);
    h_pylon = 60e-3 + 0.5 * wing_thickness;
    l_pylon = chord_root * 1.2;
    volume_spar_web = 2 * thickness_web * thickness_airfoil_relative * 2 * S_ref / (1+taper_wing)/b_wing*...
        (b_wing/2-(1-taper_wing)*b_wing/4);
    mass_spar_web = 0.93e3 * volume_spar_web;

    masses = zeros(length(parts),1);
    switch spar_diameter
        case 1e-3
            masses(strcmp(parts,'Spar Caps'),1) = 1.2e-3;
        case 1.5e-3
            masses(strcmp(parts,'Spar Caps'),1) = 2.7e-3;
        case 2e-3
            masses(strcmp(parts,'Spar Caps'),1) = 4.7e-3;
        otherwise
            masses(strcmp(parts,'Spar Caps'),1) = 1.2e-3;
    end
    masses(strcmp(parts,'Wing'),1) = 0.93e3*(0.8e-3*S_wet_wing + 0.1*V_wing);
    masses(strcmp(parts,'Tailplane'),1) = 0.93e3 * (0.4e-3*S_wet_TP + 0.1*V_TP);
    masses(strcmp(parts,'Wing Pylon'),1) = 0.93e3 * ((3.6e-3*(h_pylon - 9e-3)*(l_pylon + 20e-3)) + 22e-3 * (h_pylon + l_pylon) * 0.8e-3 + 6e-6); %kg
    masses(strcmp(parts,'Main Landing Gear Adaptor'),1) = 14e-3;
    masses(strcmp(parts,'Wheel'),1) = 6e-3;
    masses(strcmp(parts,'Motor mount'),1) = 18.5e-3;
    masses(strcmp(parts,'ESC wire'),1) = 90e-3;
    masses(strcmp(parts,'RC-Box'),1) = 250e-3;
    masses(strcmp(parts,'Wing Actuators'),1) = 23e-3;
    masses(strcmp(parts,'Fuselage'),1) = 55e-3;
    masses(strcmp(parts,'Tailplane Adaptor'),1) = 22.7e-3;%;7e-3;
    masses(strcmp(parts,'Tailplane Connector'),1) = 1e-3;
    masses(strcmp(parts,'Wing Connector'),1) = 18e-3;
    masses(strcmp(parts,'Wing Separation'),1) = 27e-3;
    masses(strcmp(parts,'Propeller'),1) = mass_propeller;
    masses(strcmp(parts,'Motor'),1) = mass_motor;
    masses(strcmp(parts,'Battery'),1) = mass_battery;
    masses(strcmp(parts,'Main Landing Gear'),1) = 18e-3; % 18g   29g   54g
    masses(strcmp(parts,'Spar web'),1) = mass_spar_web;

    mass_total = sum(masses.*quantities);
    position_cg = sum(masses.*quantities.*positions) / mass_total;
    position_cg = position_cg_target;
    break
%     if i >= 30 && (position_cg - position_cg_target) < 0
%         break
%     elseif i >= 60
%         error('Error301:noCG','Center of Gravity did not converge')
%     end
%     if abs(position_cg - position_cg_target) < 1 * mac_wing / 100 && (position_cg - position_cg_target) < 0
%         break
%     elseif position_fuselage_tip > -chord_root*0.25 - 5e-3
%         if j == 1
%             position_RC_offset_old = position_RC_offset;
%             position_cg_old = position_cg;
%             position_fuselage_tip = -chord_root*0.25 - 5e-3;
%             position_RC_offset = position_RC_offset + 5e-3;
%             
%         else
%             dposition_RC_offset = position_RC_offset - position_RC_offset_old;
%             d_position_cg = position_cg - position_cg_old;
%             position_RC_offset_old = position_RC_offset;
%             position_cg_old = position_cg;
%             position_RC_offset = position_RC_offset +...
%                 (position_cg_target - position_cg) * d_position_cg / dposition_RC_offset;
%         end
%         j = j + 1;
%     else
%         if i == 1
%             position_fuselage_old = position_fuselage_tip;
%             if position_cg - position_cg_target > 0
%                 position_fuselage_tip = position_fuselage_tip - 0.01;
%             else
%                 position_fuselage_tip = position_fuselage_tip + 0.01;
%             end
%             position_cg_old = position_cg;
%         else
%             d_position_cg = position_cg - position_cg_old;
%             d_position_fuselage_tip = position_fuselage_tip - position_fuselage_old;
%             position_fuselage_old = position_fuselage_tip;
%             position_cg_old = position_cg;
%     
%             position_fuselage_tip = position_fuselage_tip +...
%                 (position_cg_target - position_cg) * d_position_cg / d_position_fuselage_tip;
%         end
%         j = 1;
%     end
%     i = i + 1;
end
aircraftdata.Mass.mass_table.part = parts;
aircraftdata.Mass.mass_table.cg_position = positions;
aircraftdata.Mass.mass_table.quantity = quantities;
aircraftdata.Mass.mass_table.mass = masses;
aircraftdata.Mass.real_position_CG = 0.017;%position_cg;
aircraftdata.Mass.total_mass = mass_total * offsetFactor;
end