function dqdt = Dq(dtlnbeta,q,t,xc,tc,v,s,xoff,n0,Dn)
    % Find the value of dqdt along a point in the characteristic

    % interpolate to find requested value from ode45
    X = interp1(tc,xc,t);

    % get the value of dqdt
    dqdt = dtlnbeta(X,t,n0,Dn,v,s,xoff).*q;
end
