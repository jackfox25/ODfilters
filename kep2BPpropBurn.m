% ============================================================
% Function: kep2BPpropBurn
% Author: Jack Fox
% Date: 11/06/24
%     Integrator for 2BP EOM with ZOH burn prescribed
% ============================================================
% Inputs:
%        t   (1x1) time (s)
%        x   (6x1) state, same units as mu
%        u   (3x1) burn vector, frame specified in umode
%    umode   (str) input frame, can be:
%                    RIC: radial/in-track/cross-track
%                     LR: local rotating frame
%                    VEL: velocity direction (in which case u is scalar)
%                    ECI: earth-centered inertial
%       mu   (1x1) gravitational parameter of central body
% Outputs:
%     xdot   (6x1) state derivative
% ============================================================
function xdot = kep2BPpropBurn(t,x,u,umode,mu)

    r = x(1:3);
    v = x(4:6);

    kep = rv2kepler(r,v,mu);

    if strcmpi(umode,"RIC")
        % Radial / In-track / Cross-track
        u = RIC2ECI(kep.inc,kep.raan,kep.w,kep.thstar) * u;

    elseif strcmpi(umode,"LR")
        % Local rotating frame
        A = u(1:3);
        B = u(4:6);
        uLR = u(7)*(A+B*t)/norm(A+B*t);
        uRIC = [0 0 1;...
                1 0 0;...
                0 1 0] * uLR;

        u = RIC2ECI(kep.inc,kep.raan,kep.w,kep.thstar) * uRIC;

    elseif strcmpi(umode,"VEL")
        % Velocity frame
        u = u .* v/norm(v);

    elseif strcmpi(umode,"ECI")
        % u is already an ECI vector

    end

    rr = sqrt(r'*r);

    xdot = [v;
            -mu/(rr^3) * r + u];

end