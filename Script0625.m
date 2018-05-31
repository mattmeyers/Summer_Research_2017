%% Solving a basic DDE with and without AA
clear all, close all
dom = [0,5]; hdom = [-1,0]; tol = 1e-8; err = 2*tol; errVec1 = [];

init = chebfun(@(t)1,dom,'splitting','on');
hist = chebfun(@(t)1,hdom,'splitting','on');
y = join(hist,init); y0 = hist(0); k = 1;

while err > tol
    y(:,k+1) = join(hist,y0 + cumsum(chebfun(@(t)-y(t-1,k),dom,'splitting','on')));
    err = norm(y(:,k)-y(:,k+1)); errVec1(k) = err; k = k+1;
end

plot(y(:,end)), figure, semilogy(errVec1,'-o')

%% With AA
dom = [0,5]; hdom = [-1,0]; tol = 1e-8; err = 2*tol; errVec2 = [];
hist = chebfun(@(t)1,hdom,'splitting','on'); y0 = hist(0);
g = @(y) join(hist, y0 + cumsum(chebfun(@(t)-y(t-1),dom,'splitting','on')));
init = join(hist, chebfun(@(t)1,dom,'splitting','on'));

[soln,~,numIter,errVec2] = AA(g,init);

figure,plot(soln), figure, semilogy(errVec2,'-o')

%% Comparing errors
figure, semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o')
legend('Without AA', 'With AA')

%% Solving a state dependent DDE
hdom = [0,.1]; dom = [.1,5]; tol = 1e-8; err = 2*tol; errVec = [];
hist = chebfun(@(t)1./t,hdom,'splitting','on'); y0 = hist(.1);
g = @(y)join(hist,y0+cumsum(chebfun(@(t)-y(exp(1-y(t))).*y(t).^2.*exp(1-y(t)),dom,'splitting','on')));
init = join(hist,chebfun(@(t)1./(t+0.001),dom,'splitting','on'));

[soln,x,numIter,errVec]=AA(g,init);


%% Solving ODE with ODE45
fun = @(t,y) -abs(y).^2.*y+3*cos(t);
tspan = [0,20];
y0 = 0;
[t,y] = ode45(fun,tspan,y0);
plot(t,y)
