function aircraftdata = inititalComponentDrag(aircraftdata,inputdata)
%% Extract data from aircraftdata
mass = aircraftdata.Mass.total_mass;
S_ref = aircraftdata.Geometry.wing.projected_surface;
%% Extract data from inputdata
airspeed = inputdata.speed_pattern;
gravity = inputdata.gravity;
density = inputdata.density;
%% Calculation
% C_L = 2 * mass * gravity / (density * S_ref * airspeed^2);
% [C_D, ~] = dragCoefficient(aircraftdata,inputdata,airspeed,C_L);
% C_D_other = C_D / 0.35 * 0.65;
C_D_other = 0.03;
aircraftdata.Aerodynamics.drag_coefficient_other = C_D_other;
end