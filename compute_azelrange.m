function [AZ, EL, RANGE] = compute_azelrange(userECEF, satECEF, R)

    % Set up vectors with 0:Nx1
    Nsats = size(satECEF,1);
    [AZ,EL,RANGE] = deal(zeros(Nsats,1));

    % Compute user lat/lon and ECEF2ENU rotation
    userLL = ECEF2llh(userECEF,R);
    C_ECEF2ENU = ECEF2ENU(userLL(1),userLL(2));

    for i=1:Nsats
    
        % Compute relative position in ECEF
        relPos = satECEF(i,:) - userECEF;
        % Rotate relative position into ENU
        satENU = C_ECEF2ENU * relPos'; 

        % Compute AZ/EL/RANGE
        RANGE(i) = norm(relPos);
        AZ(i) = atan2d(satENU(1),satENU(2));
        EL(i) = asind(satENU(3)/RANGE(i));
    
    end

end
