clear all
close all
format long e

plot_out=0;
calc_q = 1;
disp_sanity = 0;
%% Setup q0, A: Ammplitude, sq: pulse width, qxoff: xoffset
A     =  1.0;
sq    =  2.0;
qxoff = -5.0;

q = @(x,A,s,xoff) A.*exp(-((x-xoff).^2)./s.^2);

%% Material background and RIP eta = 0 + Dn*exp(...)
s    = 5.0;
v    = 0.59;
xoff = 10.0;
n0   = 1.5;
Dn   = 0.15;

% eta
eta  = @(x,t,n0,Dn,v,s,xoff) n0 + Dn.*exp(-((x-v.*t-xoff).^2)./s.^2);

% beta = 1/eta
beta = @(x,t,n0,Dn,v,s,xoff) 1./(n0 + Dn.*exp(-((x-v.*t-xoff).^2)./s.^2));

% dt[ln(beta)]
dtlnbeta = @(x,t,n0,Dn,v,s,xoff) (2.0*Dn.*v.*(t.*v-x+xoff))./((Dn+exp(((t.*v-x+xoff).^2)./s.^2).*n0).*s.^2);

%% Grid setup
% grid 0 <= x <= 300, equivalent to PyClaw
np = 2^17;
dx = 300/np;
xf = linspace(dx/2.0,300-dx/2.0,np)';

% set intial and final times
ti = 0.0;
tf = 247.5;

% set options for ode45
options = odeset('RelTol',1e-10,'AbsTol',2.22045e-14);

%% Solve the inverse problem for the characteristics

% guess initial x
xg = xf - v.*tf;
xi = zeros(np,1);
sanity = zeros(np,1);
% scan over the values of xg
for k=1:np
    c = num2str(k);
%   trial function for fzero
    fun = @(x) dcharacteristics(beta,x,xf(k,1),ti,tf,v,s,xoff,n0,Dn,options);
    
%   find the initial point
    xi(k,1) = fzero(fun,xg(k,1));

%   get the characteristic beginning at xi
    [chars.(['t',c]),chars.(['x',c])] = ode45(@(t,x) beta(x,t,n0,Dn,v,s,xoff),[ti tf],xi(k,1),options);

%   compare if the final points are the same
    sanity(k,1) = chars.(['x',c])(end) - xf(k,1);
    if disp_sanity==1
        disp(c)
        disp(sanity(k,1))
    end

%% Solve the forward problem for q

    if calc_q==1
%       get q  initial,  q(xi,ti)
        qo = q(xi(k,1),A,sq,qxoff);

%       solve the equation along the characteristic
        [sol.(['t',c]),sol.(['q',c])] = ode45(@(t,q) Dq(dtlnbeta,q,t,chars.(['x',c]),chars.(['t',c]),v,s,xoff,n0,Dn),[ti tf],qo,options);
    end
end

%% create array with Qi and Qf
Q = zeros(np,2);
X = Q;
for i = 1:np
    c = num2str(i);
    Q(i,1) = sol.(['q',c])(1);
    Q(i,2) = sol.(['q',c])(end);
    X(i,1) = chars.(['x',c])(1);
    X(i,2) = chars.(['x',c])(end);
end
disp('Norm(1) of ode x vs given x')
disp(norm(X(:,2)-xf,1))

%% save results
basename = 'analytic';
savedir  = './results/';
save([savedir,basename,'_','centers','_all_hd_exact_',num2str(floor(np))])

%% plot
figdir   = './figures/';
if plot_out==1
    figure
    subplot(1,2,1)
    plot(chars.t,chars.x,'r--')
    xlabel('t (a.u.)')
    ylabel('x (a.u.)')
    title('Characteristic')

    subplot(1,2,2)
    plot(sol.t,sol.q,'b--')
    xlabel('t (a.u.)')
    ylabel('q (a.u.)')
    title('Solution')
end