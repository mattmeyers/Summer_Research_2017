%% Solving y''+sign(y)=0 without systems
clear all, close all
dom = [0,5]; y0 = 1; v0 = 0;
tol = 1e-8; err = 2*tol; k = 1;

t = chebfun('t', dom,'splitting','on');
y = sin(t);

while err > tol
   y(:,k+1) =  y0 + v0*t + cumsum(chebfun(@(s)(t(s)-s).*sign(y(s,k)),dom,'splitting','on'));
   err = norm(y(:,k+1)-y(:,k)); k = k+1;
end
plot(y(:,end))

%% What about with AA
clear all, close all
dom = [0,5]; y0 = 1; v0 = 0;

t = chebfun('t',dom);
init = sin(t);
g = @(y)y0 + v0*t + cumsum(chebfun(@(s)(t(s)-s).*sign(y(s)),dom,'splitting','on'));

soln = AA(g,init);