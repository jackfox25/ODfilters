function llh = ECI2llh(ECI,thetaG)

    llh = ECEF2llh(ECI2ECEF(thetaG)*ECI);

end