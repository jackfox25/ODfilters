% =========================================================================
% Function: ECEF2ECI
% Author:   Jack Fox
% Date:     02/12/25
%     DCM for transforming ECEF to ECI frame.
% =========================================================================
% Inputs:
%       thetaG   (1x1) Relative angle between ECEF and ECI frames (deg)
% Outputs:
%            C   (3x3) DCM
% =========================================================================
function C = ECEF2ECI(thetaG)
    
    C = ECI2ECEF(thetaG)';

end