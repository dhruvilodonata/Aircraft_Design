function C_L = CLforMinCD(aircraftdata,inputdata,airspeed)
mac_wing = aircraftdata.Geometry.wing.mean_aerodynamic_chord;
airfoilName_wing = aircraftdata.Geometry.wing.airfoil;

airfoil_data = inputdata.airfoil_data;
density = inputdata.density;
viscosity = inputdata.viscosity;

Re = density * airspeed * mac_wing / viscosity;

Re_list = airfoil_data.(airfoilName_wing).polar.Re_list;
polar_list = airfoil_data.(airfoilName_wing).polar.Polars;
C_L_list = zeros(1,length(Re_list));

for i = 1:length(Re_list)
    polar = polar_list{i};
    [~,index] = min(polar.CD);
    C_L_list(i) = polar.CL(index(1));
end
C_L = interp1(Re_list,C_L_list,Re);
end