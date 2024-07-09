% Calculates position of tail
% Calculates tail area and V-tail angle
% Calculates wetted area and volume of tail
% TO DO:
function [aircraftdata,error_cell] = tailSizing(aircraftdata, inputdata, error_cell)
%% Fixed Values
c_HT = 0.5;                                     % horizontal tail volume coefficient
c_VT = 0.04;                                    % vertical tail volume coefficient
l = 0.5:0.01:0.95;                               % l_HT = l_VT = l (lever arm of the tail)

%% Extract values from aircraft data struct
mac_wing = aircraftdata.Geometry.wing.mean_aerodynamic_chord;
S_ref = aircraftdata.Geometry.wing.projected_surface;
b_w = aircraftdata.Geometry.wing.span;
AR_tail = aircraftdata.Geometry.tail.aspect_ratio;
taperRatio_tail = aircraftdata.Geometry.tail.taper_ratio;

airfoil_tail = aircraftdata.Geometry.tail.airfoil;
airfoil_tail_relative_crosssection = inputdata.airfoil_data.(airfoil_tail).geo.relative_crosssection;
airfoil_tail_relative_circumference = inputdata.airfoil_data.(airfoil_tail).geo.relative_circumference;

%% Calculation
S_HT = c_HT * mac_wing * S_ref ./ l;            % area horizontal tail
S_VT = c_VT * b_w * S_ref ./ l;                 % area vertical tail
S_TP = S_VT + S_HT;                             % tail reference area

b_tail = sqrt(AR_tail*S_TP);                    % tail span
Tail_chord = S_TP/b_tail;                       % tail mean chord
thickness_TP = 0.0555 * Tail_chord;             % average thickness of tail (0.0555 per unit chord for SD8020 airfoil)

S_w_TP = 2 * S_TP;                              % wetted area tail (Could be done better)
V_TP = S_TP * thickness_TP;                     % Volume tail (Could be done better)
m_TP = 0.93e3 * (1.4e-3 * S_w_TP + 0.1 * V_TP); % mass of the tailplane
m_fuse = 55e-3 * l;                             % mass of fuselage portion behind wing
m = m_TP + m_fuse;                              % total mass that depends on tail sizing

[~, index] = min(m);                     % searches for lowest mass and also ouputs its index

l_best = l(index);
S_TP = S_TP(index);
S_HT = S_HT(index);
S_VT = S_VT(index);
b_tail = sqrt(AR_tail*S_TP);                    % tail span
mac_tail = S_TP/b_tail;                       % tail mean chord

VtailAngle = atand(S_VT / S_HT); % Defined from y-Axis

% Volume and wetted area
% How many tail sections?
chord_root_tail = 2 * mac_tail / (1 + taperRatio_tail);
tail_seperation_y = [0, b_tail/2];
tail_seperation_chord = chord_root_tail - (chord_root_tail - chord_root_tail * taperRatio_tail)/(b_tail/2)*tail_seperation_y;
tail_seperation_crosssections = airfoil_tail_relative_crosssection * tail_seperation_chord.^2;     % per tailplane
tail_seperation_crosssections_sum = 2 * sum(tail_seperation_crosssections); % per tailplane

volume_tail = airfoil_tail_relative_crosssection * (0.5 * b_tail * chord_root_tail^2 - chord_root_tail^2 * (1 - taperRatio_tail) * b_tail/2 ...
                + chord_root_tail^2 * (1 - taperRatio_tail)^2 * b_tail / 6); % Is this correct?
S_wet_tail = airfoil_tail_relative_circumference * mac_tail * b_tail + 2 * tail_seperation_crosssections_sum;

%% Write results to aircraft data struct
aircraftdata.Geometry.tail.projected_surface = S_TP;
aircraftdata.Geometry.tail.mean_aerodynamic_chord = mac_tail;
aircraftdata.Geometry.tail.span = b_tail;
aircraftdata.Geometry.tail.x_position = l_best;
aircraftdata.Geometry.tail.dihedral = VtailAngle;
aircraftdata.Geometry.tail.volume = volume_tail;
aircraftdata.Geometry.tail.wetted_surface = S_wet_tail;

% This writes Neutral Point Position and Target C.G. Position:
airspeed = 18;
C_L = 0.5;
[~,~,~,~,~,~,~,position_NP] = AVL_Call(aircraftdata,inputdata,airspeed,C_L);
aircraftdata.Aerodynamics.position_neutral_point = position_NP;
aircraftdata.Aerodynamics.target_position_CG = position_NP - 0.18 * mac_wing;
end