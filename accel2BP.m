function a = accel2BP(r,mu)

    rmag = norm(r);
    a = -mu*r/rmag^3;

end