function Q = getQ_SNC(state,params)

    Q = params.SNC.Q;

    if params.SNC.Qframe == "RIC"
        C_RIC2ECI = RIC2ECI(state);
        Q = C_RIC2ECI * Q * C_RIC2ECI';
    end

end