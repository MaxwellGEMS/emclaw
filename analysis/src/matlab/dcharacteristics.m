function y = dcharacteristics(f,x0,x1,t0,t1,v,s,xoff,n0,Dn,options)
	% solve the characteristic equation f and return the difference between
	% the calculated  and  requested values. The results is use in blacbox_*
	% to feed fzero.

	% solve characteristic equation f
    [~,x] = ode45(@(t,x) f(x,t,n0,Dn,v,s,xoff),[t0 t1],x0,options);
    
    % calculate difference
    y = x1 - x(end);
end