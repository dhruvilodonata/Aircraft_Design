function [ aircraftdata, error_cell ] = Printingtime( aircraftdata, inputdata, error_cell)
S_wet = aircraftdata.Geometry.wing.wetted_surface;
V_wing = aircraftdata.Geometry.wing.volume;
b = aircraftdata.Geometry.wing.span;
V_TP = aircraftdata.Geometry.tail.volume;
S_wet_TP = aircraftdata.Geometry.tail.wetted_surface;
b_TP = aircraftdata.Geometry.tail.span;

t_printing_wing = ((0.0008*S_wet + 0.1*V_wing)/(3.6*10^-9)) + (b/0.0003)*11;
t_printing_TP = ((0.0004*S_wet_TP + 0.1*V_TP)/(3.6*10^-9)) + (b_TP/0.0003)*11;
total_print_time = t_printing_TP + t_printing_wing;
if total_print_time > 110*3600
    error('Error302:printingTime','Printing Time exceeds 110h')
end
aircraftdata.Performance.print_time = total_print_time;