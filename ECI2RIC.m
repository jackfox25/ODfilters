function C = ECI2RIC(state)

    r = state(1:3);
    v = state(4:6);

    ihat_r = r/norm(r);
    ihat_h = cross(r,v)/norm(cross(r,v));
    ihat_th = cross(ihat_h,ihat_r);

    C = [ihat_r'; ihat_th'; ihat_h'];

end