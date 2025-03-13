% ============================================================
% Function: ECI2RIC_kep
% Author: Jack Fox
% Date: 9/20/23
%     Generates rotation matrix for transformation from 
%     inertial coordinates of central body (denoted ECI here)
%     to Radial/In-Track/Cross-Track (s/c body coordinates)
% ============================================================
% Inputs:
%      inc   inclination (deg)
%     raan   right ascension of the ascending node (deg)
%        w   argument of the periapsis (deg)
%   thstar   true anomaly (deg)
% Outputs:
%        C   direction cosine matrix 
% ============================================================
function C = ECI2RIC_kep(inc, raan, w, thstar)
    
    C_RIC2ECI = RIC2ECI_kep(inc,raan,w,thstar);
    C = C_RIC2ECI';

end