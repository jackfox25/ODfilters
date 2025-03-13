function Q_DMC = getQ_DMC(dt,state,params)

    Q = params.DMC.Q;
    B = params.DMC.B;

    if params.DMC.Qframe == "RIC"
        C_RIC2ECI = RIC2ECI(state);
        Q = C_RIC2ECI * Q * C_RIC2ECI';
    end

    B2 = B.*B;
    B3 = B2.*B;
    B4 = B3.*B;
    B5 = B4.*B;

    e1 = exp(-B*dt);
    e2 = exp(-2*B*dt);

    Qrr = Q .* ( 1./(3*B2)*dt^3 - 1./(B3)*dt^2 + 1./(B4)*dt - 2./(B4).*e1*dt + 1./(2*B5).*(1-e2) );
    Qrv = Q .* ( 1./(2*B2)*dt^2 - 1./(B3)*dt + 1./(B3).*e1*dt + 1./(B4).*(1-e1) - 1./(2*B4).*(1-e2) );
    Qrw = Q .* ( 1./(2*B3).*(1-e2) - 1./(B2).*e1*dt );
    Qvv = Q .* ( 1./(B2)*dt - 2./(B3).*(1-e1) + 1./(2*B3).*(1-e2) );
    Qvw = Q .* ( 1./(2*B2).*(1+e2) - 1./(B2).*e1 );
    Qww = Q .* ( 1./(2*B).*(1-e2) );

    Q_DMC = [Qrr Qrv Qrw;...
             Qrv Qvv Qvw;...
             Qrw Qvw Qww];

    Q_DMC(isnan(Q_DMC)) = 0;

end