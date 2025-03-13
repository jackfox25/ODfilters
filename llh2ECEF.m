% Rotates lat/lon/alt (deg) to ECEF (assuming spherical Earth)
function ECEF = llh2ECEF(llh,R)

    lat = llh(:,1);
    lon = llh(:,2);
    alt = llh(:,3);

    rho = R+alt;
    
    z = rho.*sind(lat);
    x = rho.*cosd(lat).*cosd(lon);
    y = rho.*cosd(lat).*sind(lon);
    
    ECEF = [x y z];

return