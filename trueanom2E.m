% ============================================================
% Function: trueanom2E
% Author: Jack Fox
% Date: 9/28/23
%     Converts true anomaly to eccentric anomaly.
% ============================================================
% Inputs: (units must match)
%   thstar   true anomaly (deg)
%   ecc      eccentricity
% Outputs:
%   E        eccentric anomaly (rad)
% ============================================================
function E = trueanom2E(thstar, ecc)

    E = acos( (ecc+cosd(thstar)) / (1+ecc*cosd(thstar)) );
    
    % sign check
    if mod(thstar,360) > 180
        E = -E;
        E = mod(E,2*pi);
    end



end