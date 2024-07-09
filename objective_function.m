function [aircraftdata,error_cell] = objective_function(aircraftdata,inputdata,error_cell)
    %% Precalculation
    [ aircraftdata, error_cell ] = Precalculation(aircraftdata,inputdata,error_cell);

    %% Propulsion
    [ aircraftdata, error_cell ] = Propulsion( aircraftdata, inputdata, error_cell );
    
    %% Tail Sizing
    [ aircraftdata, error_cell ] = tailSizing(aircraftdata, inputdata, error_cell);
    
    %% Preliminary Mass Calculation
    [ aircraftdata, error_cell ] = Mass( aircraftdata, inputdata, error_cell);
    
    %% Aerodynamics
    % Inititalize Component Drag
    aircraftdata = inititalComponentDrag(aircraftdata,inputdata);
    % Max possible Velocity
    [ aircraftdata, error_cell ] = maxSpeed(aircraftdata,inputdata,error_cell);
    % Lift Distribution
    [ aircraftdata, error_cell ] = liftDistribution(aircraftdata,inputdata,error_cell);
    
    %% Structure
    [ aircraftdata, error_cell ] = Structure( aircraftdata, inputdata, error_cell );
    
    %% Mass
    [ aircraftdata, error_cell ] = Mass( aircraftdata, inputdata, error_cell);
    
    %% Stall Speed Calculation
    [ aircraftdata,error_cell ] = stallSpeed( aircraftdata,inputdata,error_cell );
    
    %% Printingtime
    [ aircraftdata, error_cell ] = Printingtime( aircraftdata, inputdata, error_cell);
    
    %% Performance
    [ aircraftdata, error_cell ] = Performance( aircraftdata, inputdata, error_cell);
    
    %% Lift Distribution
    [ aircraftdata, error_cell ] = liftDistribution(aircraftdata,inputdata,error_cell);
    
    %% Points
    [ aircraftdata, error_cell ] = Points( aircraftdata, inputdata, error_cell);
end