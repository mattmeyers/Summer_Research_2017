%% Trying to solve y''+y=0 with Newton
clear all, close all

%% First with ode45 as a baseline
dom = [0,5]; IC = [1;0];
f = @(t,x) [x(2); -x(1)];
[t,x] = ode45(f,dom,IC);
plot(t,x(:,1))

%% Now with Newton
