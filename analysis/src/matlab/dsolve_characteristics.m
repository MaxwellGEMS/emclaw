% An example showing how to solve multiple characteristics for problem of interest.

close all
clear all

s = 5.0;
v = 0.61;

f = inline('1/(1.5 + 0.15*exp(-(x-v*t-10)^2/s^2))');

cc = hsv(21);
figure
hold on
for i=-10:10;
    [t,x] = ode45(@(t,x) f(s,t,v,x),[0 200],i);
    plot(t,x,'color',cc(i+11,:));
end
