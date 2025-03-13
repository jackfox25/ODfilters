% ============================================================
% Function: kepler2rv
% Author: Jack Fox
% Date: 9/20/23
%     Converts keplerian orbital elements to s/c 
%     r and v vector to.
% ============================================================
% Inputs: (units must match)
%   kepler   6x1 structure composed of fields:
%             a:  semi-major axis
%           ecc:  eccentricity
%           inc:  inclination (deg)
%          raan:  right ascension of the ascending node (deg)
%             w:  argument of the periapsis (deg)
%        thstar:  true anomaly (deg)
%   mu       GM parameter of central body
% Outputs:
%   r        3x1 central body centered inertial position vector 
%   v        3x1 central body centered inertial velocity vector
% ============================================================
function [r,v] = kepler2rv(kepler, mu)

    a = kepler.a;
    ecc = kepler.ecc;
    inc = kepler.inc;
    raan = kepler.raan;
    w = kepler.w;
    thstar = kepler.thstar;

    hmag = sqrt(mu*a*(1-ecc^2));
    rmag = hmag^2/(mu*(1+ecc*cosd(thstar)));

    r_RIC = [rmag; 0; 0];
    v_RIC = [mu/hmag*ecc*sind(thstar); mu/hmag*(1+ecc*cosd(thstar)); 0];

    r = RIC2ECI_kep(inc,raan,w,thstar) * r_RIC;
    v = RIC2ECI_kep(inc,raan,w,thstar) * v_RIC;

end