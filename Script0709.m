%% Working with the sign function
dom = [-1,1];
t = chebfun('t',dom,'splitting','on');
y1 = sign(t);
y2 = tanh(500*t);
plot([y1,y2])

%% Trying to solve y''+tanh(500y)=0 instead
clear all, close all
dom = [0,5];
g = @(x) [1 + cumsum(x(:,2)),cumsum(-tanh(500*x(:,1)))];
init = [chebfun(@(t)1,dom),chebfun(@(t)0,dom)];
[soln,x,numIter,errVec] = AA(g,init);
plot(soln)

%% Comparing to ode45
g = @(x) [1 + cumsum(x(:,2)),cumsum(-sign(100*x(:,1)))];
soln2 = AA(g,init);
