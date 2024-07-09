function inputdata = read_inputdata(fixed_parameter,variation_parameter,data_directories,airfoil_polar_naming_convention)
error_cell = cell(1,1);
%[inputdata.propulsion, error_cell] = read_propulsion_data(data_directories.propulsion, error_cell); % Propulsion (Motor, Propeller, Battery) Data
inputdata.density = 1.185;
inputdata.gravity = 9.81;
inputdata.speed_pattern = 18;
inputdata.viscosity = 17.6e-6;
inputdata.AVLOutput = false;
[inputdata.propulsion, error_cell] = import_propulsion_data(error_cell);
all_airfoils = unique([fixed_parameter.tail_airfoil, variation_parameter.wing_airfoil]); % List of all airfoils
for i_airfoil = 1:numel(all_airfoils)
    airfoil = all_airfoils{i_airfoil};
    [inputdata.airfoil_data.(airfoil).polar, error_cell] = read_airfoil_polars(airfoil, data_directories.polars, airfoil_polar_naming_convention, error_cell); % Load aerodynamic airfoil data
    [inputdata.airfoil_data.(airfoil).geo, error_cell] = airfoil_geometry(airfoil, data_directories.airfoil_geometry, fixed_parameter.spar_position, error_cell); % Load geometric airfoil data
end
end