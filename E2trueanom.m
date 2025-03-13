% ============================================================
% Function: E2trueanom
% Author: Jack Fox
% Date: 9/28/23
%     Converts eccentric anomaly to true anomaly.
% ============================================================
% Inputs: (units must match)
%   E        eccentric anomaly (rad)
%   ecc      eccentricity
% Outputs:
%   thstar   true anomaly (deg)
% ============================================================
function thstar = E2trueanom(E, ecc)

    thstar = acosd( (ecc-cos(E)) / (-1+ecc*cos(E)) );
    
    % sign check
    if mod(E,2*pi) > pi
        thstar = -thstar;
    end
    
end