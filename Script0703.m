%% Attempting Kermack-McKendrick with working AA
clear all, close all
dom = [0,1]; IC = [5,.1,1];

x1 = chebfun(5,dom);
x2 = chebfun(.1,dom);
x3 = chebfun(1,dom);
init = [x1,x2,x3];

%g1 = @(x) IC(1) + cumsum(chebfun(@(t)-x(t,1).*x(t,2),dom,'splitting','on'));
%g2 = @(x) IC(2) + cumsum(chebfun(@(t)x(t,1).*x(t,2)-x(t,2),dom,'splitting','on'));
%g3 = @(x) IC(3) + cumsum(chebfun(@(t)x(t,2),dom,'splitting','on'));
%g = @(x)[g1(x),g2(x),g3(x)];
g = @(x) [IC(1) + cumsum(-x(:,1).*x(:,2)), ...
    IC(2) + cumsum(x(:,1).*x(:,2)-x(:,2)), ...
    IC(3) + cumsum(x(:,2))];

soln = AA(g,init);

%% Let's try with ode45
dom = [0,10]; IC = [5,.1,1];
f = @(t,x) [-x(1).*x(2); x(1).*x(2)-x(2); x(2)];

[t,x] = ode45(f,dom,IC,odeset('RelTol',100,'AbsTol',100));
plot(t,x), legend('Susceptible','Infected','Recovered')

%% Using ode45 solution as initial guess for AA
clear all, close all
dom = [0,10]; IC = [5,.1,1];
f = @(t,x) [-x(1).*x(2); x(1).*x(2)-x(2); x(2)];

[t,x] = ode45(f,dom,IC,odeset('RelTol',100,'AbsTol',100));

init = [chebfun(x(:,1),dom), chebfun(x(:,2),dom), chebfun(x(:,3),dom)];
g = @(x) [IC(1) + cumsum(chebfun(@(t)-x(t,1).*x(t,2),dom,'splitting','on')), ...
    IC(2) + cumsum(chebfun(@(t)x(t,1).*x(t,2)-x(t,2),dom,'splitting','on')), ...
    IC(3) + cumsum(chebfun(@(t)x(t,2),dom,'splitting','on'))];

soln = AA(g,init);
plot(soln)