% ============================================================
% Function: soloInputProp
% Author: Jack Fox
% Date: 12/08/24
%     General integrator for [r v]' state with prescribed 
%     ZOH input vector.
% ============================================================
% Inputs:
%        t   (1x1) time (s)
%        x   (6x1) state
%        u   (3x1) burn vector
% Outputs:
%     xdot   (6x1) state derivative
% ============================================================
function xdot = soloInputProp(t,x,u)

    v = x(4:6);

    xdot = [v;
            u];

end