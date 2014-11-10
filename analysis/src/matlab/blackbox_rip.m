clear all
close all
format long e
disp('running')
plot_out=0;
calc_q = 1;
disp_sanity = 1;
%% Setup q0, A: Ammplitude, sq: pulse width, qxoff: xoffset
A     =  1.0;
sq    =  2.0;
qxoff = 10.0;

q = @(x,A,s,xoff,dx) A.*exp(-((x-xoff).^2)./s.^2);
% q = @(x,A,s,xoff,dx) A.*sqrt(pi)*s*(erf(((dx/2.0) + (xoff - x))/s)+erf(((dx/2.0) - (xoff - x))/s))/(2.0*dx);

%% Material background and RIP eta = 0 + Dn*exp(...)
s    = 5.0;
v    = 0.59;
xoff = 25.0;
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
np = 2^16;
Lx = 300.0;
dx = Lx/np;
mp = 2*np+1;
xf = linspace(0,Lx,mp)';

% set intial and final times
ti = 0.0;
tf = 274.5;

% set options for ode45
options = odeset('RelTol',1e-10,'AbsTol',2.22045e-14);

%% Solve the inverse problem for the characteristics

% guess initial x
xg = xf - v.*tf;
xi = zeros(np,1);
sanity = zeros(np,1);
% scan over the values of xg
for k=1:mp
    c = num2str(k);
%   trial function for fzero
    fun = @(x) dcharacteristics(beta,x,xf(k,1),ti,tf,v,s,xoff,n0,Dn,options);
    
%   find the initial point
    xi(k,1) = fzero(fun,xg(k,1));

%   get the characteristic beginning at xi
    [chars.(['t',c]),chars.(['x',c])] = ode45(@(t,x) beta(x,t,n0,Dn,v,s,xoff),[ti tf],xi(k,1),options);

%   compare if the final points are the same
    sanity(k,1) = chars.(['x',c])(end) - xf(k,1);
    if (disp_sanity==1 && mod(k,100)==0)
        disp(c)
        disp(sanity(k,1))
    end

%% Solve the forward problem for q

    if calc_q==1
%       get q  initial,  q(xi,ti)
        qo = q(xi(k,1),A,sq,qxoff,dx);

%       solve the equation along the characteristic
        [sol.(['t',c]),sol.(['q',c])] = ode45(@(t,q) Dq(dtlnbeta,q,t,chars.(['x',c]),chars.(['t',c]),v,s,xoff,n0,Dn),[ti tf],qo,options);
    end
end

%% create array with Qi and Qf
Q = zeros(mp,2);
X = Q;
for i = 1:mp
    c = num2str(i);
    Q(i,1) = sol.(['q',c])(1);
    Q(i,2) = sol.(['q',c])(end);
    X(i,1) = chars.(['x',c])(1);
    X(i,2) = chars.(['x',c])(end);
end

qs(:,1) = (1/6.0)*(Q(1:2:(end-2),1)+4.0*Q(2:2:(end-1),1)+Q(3:2:end,1));
qs(:,2) = (1/6.0)*(Q(1:2:(end-2),2)+4.0*Q(2:2:(end-1),2)+Q(3:2:end,2));
xs(:,1) = (1/6.0)*(X(1:2:(end-2),1)+4.0*X(2:2:(end-1),1)+X(3:2:end,1));
xs(:,2) = (1/6.0)*(X(1:2:(end-2),2)+4.0*X(2:2:(end-1),2)+X(3:2:end,2));

disp('Norm(1) of ode x vs given x')
disp(norm(X(:,2)-xf,1))

%% save results
basename = '_rip';
savedir  = './results/';
save([savedir,basename,'_nc_',num2str(floor(np))])

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