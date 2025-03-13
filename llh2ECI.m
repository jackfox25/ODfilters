function ECI = llh2ECI(llh,thetaG,R)
    
    ECI = ECEF2ECI(thetaG) * llh2ECEF(llh,R)';

end