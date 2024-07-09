% Calculates the maximum 3D lift coefficient from given airspeed and wing
% geometry
% Inputs: Airspeed (can be an array)
% Outputs: Maximum 3D-Lift Coefficient (same length as Airspeed)

function C_L_max = max3DLiftCoefficient(aircraftdata, inputdata, airspeedArray)
C_L_max = zeros(1,length(airspeedArray));
for i = 1:length(airspeedArray)
    airspeed = airspeedArray(i);
    C_L = 1;
    [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,airspeed, C_L);
    if stall
        while true
            C_L = C_L - 0.1;
            [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,airspeed, C_L);
            if ~stall
                break
            elseif C_L <= 0
                error('Error103:CLmax','No valid C_L without stall found')
            end
        end
        while true
            C_L = C_L + 0.02;
            [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,airspeed, C_L);
            if stall
                C_L_max(i) = C_L - 0.02;
                break
            end
        end
    else
        while true
            C_L = C_L + 0.1;
            [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,airspeed, C_L);
            if stall
                break
            elseif C_L > 2
                error('Error103:CLmax','No stall below C_L = 2')
            end
        end
        while true
            C_L = C_L - 0.02;
            [~,~,~,~,~,~,stall] = AVL_Call(aircraftdata,inputdata,airspeed, C_L);
            if ~stall
                C_L_max(i) = C_L;
                break
            end
        end
    end
end