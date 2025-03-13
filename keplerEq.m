function [E, iter] = keplerEq (M0, M, e, tol)
%%% Newton-Raphson Method for Elliptical orbits
g  =  @(E) E-e*sin(E) - M;  % Function to find the root of
dg =  @(E) 1-e*cos(E);      % Derivative of the function
 
E = M0;                  % Initial guess
err = g(E);              % Error of guess

iter = 0;
MAX_ITER = 10e5;

while (abs(err)>tol && iter < MAX_ITER)
    E = E - g(E)/dg(E);
    err = g(E);
    iter = iter + 1;
end


