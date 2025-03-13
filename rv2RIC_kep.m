function [rRIC, vRIC] = rv2RIC_kep(rECI,vECI,mu)

    k = rv2kepler(rECI,vECI,mu);

    C_ECI2RIC = ECI2RIC(k(3),k(4),k(5),k(6));
    
    rRIC = C_ECI2RIC * rECI;
    vRIC = C_ECI2RIC * vECI;

end