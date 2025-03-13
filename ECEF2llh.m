% Rotates ECEF to lat/lon/alt (deg) (assuming spherical Earth)
function llh = ECEF2llh(ECEF,R)

    x = ECEF(1);
    y = ECEF(2);
    z = ECEF(3);

    rho = norm(ECEF);

    alt = rho-R;
    lat = asind(z/rho);
    lon = atan2d(y,x);
    
    if lon > 180
        lon = lon-360;
    elseif lon < -180
        lon = lon+360;
    end


    llh = [lat;lon;alt];

end