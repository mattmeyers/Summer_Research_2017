%% Trying some Newton iterations: y'' + y = 0
clear all, close all
dom = [0,5];
N = chebop(@(t,u,v)[u-volt(@(t,s)1,v);v+volt(@(t,s)1,u)],dom);
yh1 = chebfun(@(t)1,dom); yh2 = chebfun(@(t)1,dom);
r = [1+cumsum(yh2)-yh1;1-cumsum(yh1)-yh2];

u = N\r;

u = [u{1}+chebfun(@(t)1,dom),u{2}+chebfun(@(t)1,dom)];
plot(u)

%% The nonlinear equation: y''+y^2=0
clear all, close all
dom = [0,2]; tol = 1e-8; err = 2*tol; k = 1;
y0 = 1; v0 = 0; uh = chebfun(@(t)cos(t),dom); vh = chebfun(@(t)-sin(t),dom);
N = chebop(@(t,d,e)[d-volt(@(t,s)1,e);e+2*volt(@(t,s)uh(t),d)],dom);
r = [y0+cumsum(vh)-uh;v0-cumsum(uh.^2)-vh];

x = N\r;
x = [uh+x{1}, vh+x{2}];

while err > tol
    uh = x(:,1); vh = x(:,2);
    N = chebop(@(t,d,e)[d-volt(@(t,s)1,e);e+volt(@(t,s)2*uh(t),d)],dom);
    r = [y0+cumsum(vh)-uh;v0+cumsum(-uh.^2)-vh];

    xn = N\r;
    xn = [uh+xn{1}, vh+xn{2}];
    err = max(norm(x(:,1)-xn(:,1)),norm(x(:,2)-xn(:,2))); x = xn; k = k+1;
end

%% The nonlinear equation: y'+y^2=0
clear all, close all
dom = [0,2]; tol = 1e-8; err = 2*tol; k = 1;
y0 = 1; yh = chebfun(@(t)exp(-t),dom);
N = chebop(@(t,d)d+2*volt(@(t,s)exp(-t),d),dom);
r = y0-cumsum(yh.^2)-yh;

y = N\r;
y = yh+y;

while err > tol
    yh = y(:,1);
    N = chebop(@(t,d)d+2*volt(@(t,s)yh(t),d),dom);
    r = y0-cumsum(yh.^2)-yh;

    yn = N\r;
    yn = yh+yn;
    err = norm(y-yn); y = yn; k = k+1;
end