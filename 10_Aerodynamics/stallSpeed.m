% Calculates stall speed from given aircraft data
% TO DO:
% 
function [aircraftdata,error_cell] = stallSpeed(aircraftdata,inputdata,error_cell)
%% Extract data from inputdata
density = inputdata.density;         %density at which we operate
g = inputdata.gravity;

%% Extract data from aircraft data struct
S_ref = aircraftdata.Geometry.wing.projected_surface;
m = aircraftdata.Mass.total_mass;

%% Calculation
V = 9;
q = 0.5 * density * V^2;
C_L = m * g / (q * S_ref);
if C_L <= 2
    [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,V,C_L);
else
    stall = true;
end
if stall
    V = V + 1;
    while true
        q = 0.5 * density * V^2;
        C_L = m * g / (q * S_ref);
        if C_L > 2
            V = V + 1;
            continue
        end
        [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,V,C_L);

        if ~stall
            V_stall = V;
            break
        elseif V > 10
            error('Error104:stallSpeed','Stall speed is at least 10 m/s')
        else
            V = V + 1;
        end
    end

    % For higher fidelity
    for i = 1:10
        V = V_stall - i * 0.1;
        q = 0.5 * density * V^2;
        C_L = m * g / (q * S_ref);
        if C_L > 2
            V_stall = V + 0.1;
            break
        end
        [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,V,C_L);
        if stall
            V_stall = V + 0.1;
            break
        end
    end
else
    V = V - 1;
    while true
        q = 0.5 * density * V^2;
        C_L = m * g / (q * S_ref);
        if C_L > 2
            V_stall = V + 1;
            break
        end
        [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,V,C_L);
    
        if ~stall
            V = V - 1;
        else
            V_stall = V + 1;
            break
        end
    end
    % For higher fidelity
    for i = 1:10
        V = V_stall - i * 0.1;
        q = 0.5 * density * V^2;
        C_L = m * g / (q * S_ref);
        if C_L > 2
            V_stall = V + 0.1;
            break
        end
        [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,V,C_L);
        if stall
            V_stall = V + 0.1;
            break
        end
    end
end

aircraftdata.Performance.stall_speed = V_stall;
end

