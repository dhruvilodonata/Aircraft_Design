function [T, Q, RPM, P_tot,P_shaft,P_motor, eta_tot, eta_prop, eta_motor, eta_ESC, I_motor, V_motor, I_ESC, V_ESC] = interpolate_max_throttle_from_lookup(propulsion_lookup,ID_P,ID_M,ID_B,airspeed)

    % If airspeed is in between two integer values
    if airspeed > 0 && ~(mod(airspeed,1) == 0)
        
        % Get integer values of the airspeed next to the given airspeed
        airspeed1 = airspeed - mod(airspeed,1);
        airspeed2 = airspeed - mod(airspeed,1) + 1; 

        % Extract both airspeed tables from interpolation structures
        lookup1 = propulsion_lookup.prop(ID_P).motor(ID_M).battery(ID_B).airspeed(airspeed1 + 1).lookup_table;
        lookup2 = propulsion_lookup.prop(ID_P).motor(ID_M).battery(ID_B).airspeed(airspeed2 + 1).lookup_table;
        
        % Set airspeed interpolation points
        airspeed_pts = [airspeed1 airspeed2];
        
        % First interpolate the wanted value in both airspeed tables with
        % the given thrust. Then interpolate those values with the given
        % airspeed.

        T1 = lookup1.T(end); T2 = lookup2.T(end);
        T = interp1(airspeed_pts, [T1 T2], airspeed);

        Q1 = lookup1.Q(end); Q2 = lookup2.Q(end);
        Q = interp1(airspeed_pts, [Q1 Q2], airspeed);

        RPM1 = lookup1.RPM(end); RPM2 = lookup2.RPM(end);
        RPM = interp1(airspeed_pts, [RPM1 RPM2], airspeed);

        P_tot1 = lookup1.P_tot(end); P_tot2 = lookup2.P_tot(end);
        P_tot = interp1(airspeed_pts, [P_tot1 P_tot2], airspeed);

        P_shaft1 = lookup1.P_shaft(end); P_shaft2 = lookup2.P_shaft(end);
        P_shaft = interp1(airspeed_pts,[P_shaft1 P_shaft2],airspeed);

        P_motor1 = lookup1.P_motor(end); P_motor2 = lookup2.P_motor(end);
        P_motor = interp1(airspeed_pts,[P_motor1 P_motor2],airspeed);

        eta_tot1 = lookup1.eta_tot(end); eta_tot2 = lookup2.eta_tot(end);
        eta_tot = interp1(airspeed_pts,[eta_tot1 eta_tot2],airspeed);

        eta_prop1 = lookup1.eta_prop(end); eta_prop2 = lookup2.eta_prop(end);
        eta_prop = interp1(airspeed_pts,[eta_prop1 eta_prop2],airspeed);

        eta_motor1 = lookup1.eta_motor(end); eta_motor2 = lookup2.eta_motor(end);
        eta_motor = interp1(airspeed_pts,[eta_motor1 eta_motor2],airspeed);

        eta_ESC1 = lookup1.eta_ESC(end); eta_ESC2 = lookup2.eta_ESC(end);
        eta_ESC = interp1(airspeed_pts,[eta_ESC1 eta_ESC2],airspeed);

        I_motor1 = lookup1.I_motor(end); I_motor2 = lookup2.I_motor(end);
        I_motor = interp1(airspeed_pts,[I_motor1 I_motor2],airspeed);

        V_motor1 = lookup1.V_motor(end); V_motor2 = lookup2.V_motor(end);
        V_motor = interp1(airspeed_pts,[V_motor1 V_motor2],airspeed);

        I_ESC1 = lookup1.I_ESC(end); I_ESC2 = lookup2.I_ESC(end);
        I_ESC = interp1(airspeed_pts,[I_ESC1 I_ESC2],airspeed);

        V_ESC1 = lookup1.V_ESC(end); V_ESC2 = lookup2.V_ESC(end);
        V_ESC = interp1(airspeed_pts,[V_ESC1 V_ESC2],airspeed);

    % If the given airspeed is already an integer just interpolate all values from the corresponding lookup-table    
    elseif mod(airspeed,1) == 0
        lookup = propulsion_lookup.prop(ID_P).motor(ID_M).battery(ID_B).airspeed(airspeed + 1).lookup_table;
        T = lookup.T(end);

        Q = lookup.Q(end);

        RPM = lookup.RPM(end);

        P_tot = lookup.P_tot(end);

        P_shaft = lookup.P_shaft(end);

        P_motor = lookup.P_motor(end);

        eta_tot = lookup.eta_tot(end);

        eta_prop = lookup.eta_prop(end);

        eta_motor = lookup.eta_motor(end);

        eta_ESC = lookup.eta_ESC(end);

        I_motor = lookup.I_motor(end);

        V_motor = lookup.V_motor(end);

        I_ESC = lookup.I_ESC(end);

        V_ESC = lookup.V_ESC(end);
    end

end

