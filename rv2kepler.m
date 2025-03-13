% ============================================================
% Function: rv2kepler
% Author: Jack Fox
% Date: 9/15/23
%     Converts s/c r and v vector to keplerian
%     orbital elements.
% ============================================================
% Inputs: (units must match)
%   r        3x1 central body centered inertial position vector 
%   v        3x1 central body centered inertial velocity vector
%   mu       GM parameter of central body
% Outputs:
%   kepler   structure with fields:
%             a:  semi-major axis
%           ecc:  eccentricity
%           inc:  inclination (deg)
%          raan:  right ascension of the ascending node (deg)
%             w:  argument of the periapsis (deg)
%        thstar:  true anomaly (deg)
%   hmag     angular momentum magnitude
%   sme      specific mechanical energy
% ============================================================
function [kepler, hmag, sme, evnorm] = rv2kepler(r,v,mu)

    rmag = norm(r);
    vmag = norm(v);
    
    sme = 0.5*vmag^2-mu/rmag; % specific energy
    a = -mu/(2*sme);
    
    h = cross(r,v); % specific angular momentum
    hmag = norm(h);
    hhat = h/hmag;
    
    eccvec = cross(v,h)/mu-r/rmag; % eccentricity vector
    ecc = norm(eccvec); % eccentricity
    evnorm = eccvec./ecc;

    inc = acosd(hhat(3)) + eps; % inclination [deg]
    
    n = cross([0;0;1],h); % line of nodes
    raan = acosd(dot([1;0;0],n)/norm(n)) * n(2)/abs(n(2));
    
    if inc == 0 && ecc ~= 0
        w = acosd(dot([1;0;0],eccvec)/norm(eccvec));
        if dot(eccvec,[0;1;0]) < 0
            w = -w;
        end
    else
        w = acosd(dot(n,eccvec)/(ecc*norm(n))) * eccvec(3)/abs(eccvec(3));
    end


    thstar = acosd(dot(r,eccvec)/(rmag*ecc)) * dot(r,v)/abs(dot(r,v));
    
    kepler.a = a;
    kepler.ecc = ecc;
    kepler.inc = inc;
    kepler.raan = raan;
    kepler.w = w;
    kepler.thstar = thstar;

end