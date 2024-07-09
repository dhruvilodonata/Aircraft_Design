% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                             Cd_Cm_of_Cl.m
%
%      L E H R S T U H L   F U E R   L U F T F A H R T S Y S T E M E
%        I N S T I T U T E   O F   A I R C R A F T   D E S I G N
%                      Web: www.lls.lrg.tum.de
%                        ______
%                          |  |    |\  /|
%                          |  |    | \/ |
%                          |  |____|    |
%                Technische Universitaet Muenchen TUM
%
% (c) 2019 by Institute of Aircraft Design of TUM
%                       All Rights Reserved
%
%       P R A K T I K U M   F L U G Z E U G E N T W U R F   1 9 / 20
%                             __   __   __   __ 
%                        |   |__| |__| |__| |__|
%                        |   |  | |\   |    |  |
%                        |__ |  | | \  |    |  |
%                        
%                 L L S   Advanced Research Project Agency
%
% Enter:        -
% Module:       -
% Version:      1.0
% -------------------------------------------------------------------------
% DESCRIPTION:
% Interpolates the drag and pitching moment coefficients for an airfoil at
% a certain lift coefficient and reynolds number
% -------------------------------------------------------------------------
% DATE          VERSION     AUTHOR              CHANGELOG
% 15/10/2019    v1.0        Thomas Seren        initial release
% -------------------------------------------------------------------------
% Input variables:
% NAME          UNIT    DESCRIPTION
% Airfoil_Data_Struct  -   Airfoil Struct from read_airfoil_polars.m
% Cl     -   lift coefficient
% Reynolds_number - Reynolds number
% -------------------------------------------------------------------------
% Output variables:
% NAME          UNIT    DESCRIPTION
% Cd   -  Drag coefficient
% Cm   - Pitching moment coefficient
% -------------------------------------------------------------------------
% Limitations:  -
% Assumptions:  -
% -------------------------------------------------------------------------
% Example:      [Cd, Cm] = Cd_Cm_of_Cl(Airfoil_Data_Struct, .5, 123e3)
% -------------------------------------------------------------------------
% Called functions:
% interp1qr.m from matlab fileexchange https://www.mathworks.com/matlabcentral/fileexchange/43325-quicker-1d-linear-interpolation-interp1qr
% -------------------------------------------------------------------------
% Error Management
% ERROR MESSAGE                            DESCRIPTION
% -                                         -
% -------------------------------------------------------------------------
% KNOWN ISSUES:
% - 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [Cd, Cm] = Cd_Cm_of_Cl(Airfoil_Data_Struct, Cl, Reynolds_number)
    [~,i_Re] = histc(Reynolds_number, Airfoil_Data_Struct.Re_list); % find the next lower Reynlods number that is available
    Cds = NaN(2,1);
    Cms = NaN(2,1);
    Res = Airfoil_Data_Struct.Re_list(i_Re+[0,1]).'; %create new Reynlods number array for interpolation
    for ii = 1:2
        cur_Polar = Airfoil_Data_Struct.Polars{ii-1+i_Re}; %select Polar
        Cl_vector = cur_Polar.CL; % get CL_vector
        Cd_vector = cur_Polar.CD; % get Cd_vector
        Cm_vector = cur_Polar.Cm; % get Cm_vector
        Cds(ii) = interp1qr(Cl_vector,Cd_vector,Cl); % interpolate Cd for one Cl
        Cms(ii) = interp1qr(Cl_vector,Cm_vector,Cl); % interpolate Cm for one Cl
    end
    Cd = interp1qr(Res,Cds,Reynolds_number); % interpolate Cd for the required Reynolds number
    Cm = interp1qr(Res,Cms,Reynolds_number); % interpolate Cm for the required Reynolds number
end

