clear all
close all
format long e

plot_out=1;
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
np = 2^14;
dx = 300/np;
x  = linspace(dx/2.0,300-dx/2.0,np);

% select point of interest, xf
xf = x(8465);

% set intial and final times
ti = 0.0;
tf = 247.5;

% set options for ode45
options = odeset('RelTol',2.22045e-14,'AbsTol',5e-16);

%% Solve the inverse problem for the characteristics
% guess initial x
xg = xf - v.*tf;

% trial function for fzero
fun = @(x) dcharacteristics(beta,x,xf,ti,tf,v,s,xoff,n0,Dn,options);

% find the initial point
xi = fzero(fun,xg);

% get the characteristic beginning at xi
[chars.t,chars.x] = ode45(@(t,x) beta(x,t,n0,Dn,v,s,xoff),[ti tf],xi,options);

% compare if the final points are the same
sanity = chars.x(end) - xf;
disp(sanity)

%% Solve the forward problem for q
% get q  initial,  q(xi,ti)
qo = q(xi,A,sq,qxoff);

% solve the equation along the characteristic
[sol.t,sol.q] = ode45(@(t,q) Dq(dtlnbeta,q,t,chars.x,chars.t,v,s,xoff,n0,Dn),[ti tf],qo,options);

%% plot
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