function xdot = kep2BPprop(t,x,mu)

    r = x(1:3);
    v = x(4:6);

    rr = sqrt(r'*r);

    xdot = [v;
            -mu/(rr^3) * r];

end