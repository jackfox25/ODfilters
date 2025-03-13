function C_ECEF2ENU = ECEF2ENU(lat_deg,lon_deg)

    slon = sind(lon_deg);
    clon = cosd(lon_deg);

    slat = sind(lat_deg);
    clat = cosd(lat_deg);

    C_ECEF2ENU = [-slon       clon      0   ;...
                  -slat*clon -slat*slon clat;...
                   clat*clon  clat*slon slat];



end