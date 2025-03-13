% =========================================================================
% Function: ECI2ECEF
% Author:   Jack Fox
% Date:     02/12/25
%     DCM for transforming ECI to ECEF frame.
% =========================================================================
% Inputs:
%       thetaG   (1x1) Relative angle between ECEF and ECI frames (deg)
% Outputs:
%            C   (3x3) DCM
% =========================================================================
function C = ECI2ECEF(thetaG)

    C = [ cosd(thetaG) sind(thetaG) 0;...
         -sind(thetaG) cosd(thetaG) 0;...
          0           0           1];

end