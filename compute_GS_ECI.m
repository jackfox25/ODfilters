function GS_ECI = compute_GS_ECI(GS_llh,n_GS,R,thetaG)

    GS_ECI = zeros(n_GS,6);
    for i=1:n_GS

        lat = GS_llh(i,1);
        lon = GS_llh(i,2);

        GS_ECI(i,1:3) = llh2ECI(GS_llh(i,:),thetaG,R);

        % Compute inertial velocity at the given latitude
        vmag = 2*pi*R/86400*cosd(lat);

        % Rotate to correct vector direction
        GS_ECI(i,4:6) = ECEF2ECI(lon+thetaG) * [0; vmag; 0];

    end
    

end