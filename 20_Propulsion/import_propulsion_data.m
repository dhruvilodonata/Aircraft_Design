function [data, error_cell] = import_propulsion_data(error_cell)
propulsion_lookup = load("80_Databases/Propulsion/propulsion_lookup.mat");
data.propulsion_lookup = propulsion_lookup.propulsion_lookup;
Propulsion_Data_Struct = load("80_Databases/Propulsion/Propulsion_Data_Struct.mat");
data.Propulsion_Data_Struct = Propulsion_Data_Struct.Propulsion_Data_Struct;
end