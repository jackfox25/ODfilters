% ============================================================
% Function: RIC2ECI_kep
% Author: Jack Fox
% Date: 9/20/23
%     Generates rotation matrix for transformation from 
%     Radial/In-Track/Cross-Track (s/c body coordinates)
%     to inertial coordinates of central body (denoted ECI here)
% ============================================================
% Inputs:
%      inc   inclination (deg)
%     raan   right ascension of the ascending node (deg)
%        w   argument of the periapsis (deg)
%   thstar   true anomaly (deg)
% Outputs:
%        C   direction cosine matrix 
% ============================================================
function C = RIC2ECI_kep(inc, raan, w, thstar)
    
    if isnan(inc)
        inc = 0;
    end
    if isnan(raan)
        raan = 0;
    end
    if isnan(w)
        w = 0;
    end
    if isnan(thstar)
        thstar = 0;
    end

    th = w + thstar;

    cr = cosd(raan); sr = sind(raan);
    ci = cosd(inc); si = sind(inc);
    ct = cosd(th); st = sind(th);

    C = [cr*ct-sr*ci*st   -cr*st-sr*ci*ct    sr*si;...
         sr*ct+cr*ci*st   -sr*st+cr*ci*ct   -cr*si;...
         si*st             si*ct             ci];

end