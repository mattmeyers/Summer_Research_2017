%% EX1; y' + abs(y)*y=0 PI - Doesn't work
clear all, close all
tic
tol = 1e-8; err = 2*tol; k = 1; errVec1 = [];
dom = [0,5]; y0 = 1; y = chebfun(@(t)1,dom,'splitting','on');

y(:,2) = y0 - cumsum(abs(y).*y); 
err = norm(y(:,1)-y(:,2)); errVec1 = err; k = k+1;

while err > tol
   y(:,k+1) = y0 - cumsum(abs(y(:,k)).*y(:,k)); 
   err = norm(y(:,k)-y(:,k+1)); errVec1(k) = err; k = k+1;
end
t1 = toc;

%% PIAA - Does work
tic
g = @(y) y0 + cumsum(-abs(y).*y);
init = chebfun(@(t)exp(-t),dom);
[soln,~,numIter,errVec2] = AA(g,init);
t2 = toc;

%% Newton
tic
dom = [0,5]; tol = 1e-8; err = 2*tol; k = 1; errVec3 = [];
y0 = 1; yh = chebfun(@(t)1,dom);
N = chebop(@(t,d)d+volt(@(t,s)abs(yh(t)),d)+volt(@(t,s)yh(t).*sign(yh(t)),d),dom);
r = y0-cumsum(abs(yh).*yh)-yh;

y = N\r;
y = yh+y;
errVec3 = norm(yh-y);

while err > tol
    yh = y(:,1);
    N = chebop(@(t,d)d+volt(@(t,s)abs(yh(t)),d)+volt(@(t,s)yh(t).*sign(yh(t)),d),dom);
    r = y0-cumsum(abs(yh).*yh)-yh;

    yn = N\r;
    yn = yh+yn;
    err = norm(y-yn); errVec3(k+1) = err; y = yn; k = k+1;
end
t3 = toc;

%% Comparing errors
semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o'),semilogy(errVec3,'-o')
title('y''+abs(y)*y=0, y(0)=1: Error'),xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA','Newton')

%% With Chebfun
tic
N = chebop(@(t,y)diff(y)+abs(y).*y,dom);
y = N\0;
t4 = toc;

%% With ode45
tic
f = @(t,y) -abs(y).*y;
[t,y] = ode45(f,dom,y0);
t5 = toc
    

%% Ex2: y'+y^2=0 PI
clear all, close all
tol = 1e-8; err = 2*tol; k = 1; errVec1 = [];
dom = [0,2]; y0 = 1; y = chebfun(@(t)1,dom,'splitting','on');

y(:,2) = y0 - cumsum(y.^2); 
err = norm(y(:,1)-y(:,2)); errVec1 = err; k = k+1;

while err > tol
   y(:,k+1) = y0 - cumsum(y(:,k).^2); 
   err = norm(y(:,k)-y(:,k+1)); errVec1(k) = err; k = k+1;
end
plot(y(:,end))

%% PIAA
g = @(y) y0 + cumsum(-y.^2);
init = chebfun(@(t)1,dom);
[soln,~,numIter,errVec2] = AA(g,init);

%% Newton
dom = [0,2]; tol = 1e-8; err = 2*tol; k = 1; errVec3 = [];
y0 = 1; yh = chebfun(@(t)1,dom);
N = chebop(@(t,d)d+2*volt(@(t,s)yh(t),d),dom);
r = y0-cumsum(yh.^2)-yh;

y = N\r;
y = yh+y;
errVec3 = norm(yh-y);

while err > tol
    yh = y(:,1);
    N = chebop(@(t,d)d+2*volt(@(t,s)yh(t),d),dom);
    r = y0-cumsum(yh.^2)-yh;

    yn = N\r;
    yn = yh+yn;
    err = norm(y-yn); errVec3(k+1) = err; y = yn; k = k+1;
end

%% Comparing Error
semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o'), semilogy(errVec3,'-o')
title('y''+y^2=0, y(0)=1: Error'),xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA','Newton')


%% Ex3: y' = -(1/5)y(t-1) PI
clear all, close all
hdom = [-1,0]; dom = [0,4]; y0 = 1;
yhist = chebfun(@(t)(-1).^floor(-4*t),hdom,'splitting','on');
yinit = chebfun(@(t) 1,dom);
y = join(yhist,yinit);

tol = 1e-8; err = 1; k = 1; errVec1 = [];

y(:,2) = join(yhist,y0-(1/5)*cumsum(chebfun(@(s)y(s-1),dom,'splitting','on')));
err = norm(y(:,2)-y(:,1)); errVec1(1) = err; k = k+1;

while err > tol
    y(:,k+1) = join(yhist,y0 - (1/5)*cumsum(chebfun(@(s)y(s-1,k),dom,'splitting','on')));
    err = norm(y(:,k+1)-y(:,k));
    errVec1(k) = err;
    k = k+1;
end

%% PIAA
g = @(y) join(yhist,y0-(1/5)*cumsum(chebfun(@(s)y(s-1),dom,'splitting','on')));
init = join(yhist,yinit);
[soln,~,numIter,errVec2] = AA(g,init);

%% Comparing errors
semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o')
title('y''+(1/5)*y(t-1)=0, \phi(t)=(-1)^{floor(4t)}: Error'),xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA')


%% Ex4: y'(t)+y(t-1)=0, phi=1 PI
clear all, close all
hdom = [-1,0]; dom = [0,4]; y0 = 1;
yhist = chebfun(@(t)1,hdom,'splitting','on');
yinit = chebfun(@(t) 1,dom);
y = join(yhist,yinit);

tol = 1e-8; err = 1; k = 1; errVec1 = [];

y(:,2) = join(yhist,y0-cumsum(chebfun(@(s)y(s-1),dom,'splitting','on')));
err = norm(y(:,2)-y(:,1)); errVec1(1) = err; k = k+1;

while err > tol
    y(:,k+1) = join(yhist,y0 - cumsum(chebfun(@(s)y(s-1,k),dom,'splitting','on')));
    err = norm(y(:,k+1)-y(:,k));
    errVec1(k) = err;
    k = k+1;
end

%% PIAA
g = @(y) join(yhist,y0-cumsum(chebfun(@(s)y(s-1),dom,'splitting','on')));
init = join(yhist,yinit);
[soln,~,numIter,errVec2] = AA(g,init);

%% Comparing errors
semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o')
title('y''+y(t-1)=0, \phi(t)=1: Error'),xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA')


%% Ex5: y'-sin(t^2)y=0 PI
clear all, close all
tol = 1e-8; err = 2*tol; k = 1; errVec1 = [];
dom = [0,8]; y0 = 1; y = chebfun(@(t)1,dom,'splitting','on');
t = chebfun('t',dom);

y(:,2) = y0 + cumsum(sin(t.^2.).*y); 
err = norm(y(:,1)-y(:,2)); errVec1 = err; k = k+1;

while err > tol
   y(:,k+1) = y0 + cumsum(sin(t.^2).*y(:,k)); 
   err = norm(y(:,k)-y(:,k+1)); errVec1(k) = err; k = k+1;
end
plot(y(:,end))

%% PIAA
g = @(y) y0 + cumsum(sin(t.^2).*y);
init = chebfun(@(t)1,dom);
[soln,~,numIter,errVec2] = AA(g,init);

%% Newton
dom = [0,8]; tol = 1e-8; err = 2*tol; k = 1; errVec3 = [];
y0 = 1; yh = chebfun(@(t)1,dom);
N = chebop(@(t,d)d-volt(@(t,s)sin(s.^2),d),dom);
r = y0+cumsum(sin(t.^2).*yh)-yh;

y = N\r;
y = yh+y;
errVec3 = norm(yh-y);
%%
while err > tol
    yh = y(:,1);
    N = chebop(@(t,d)d-volt(@(t,s)sin(s.^2),d),dom);
    r = y0+cumsum(sin(t.^2).*yh)-yh;

    yn = N\r;
    yn = yh+yn;
    err = norm(y-yn); errVec3(k+1) = err; y = yn; k = k+1;
end

%% Comparing Error
semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o'), semilogy(errVec3,'-o')
title('y''-sin(t^2)y=0, y(0)=1: Error'),xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA','Newton')


%% Returning to Kermack McKendrick
clear all, close all
dom = [0,1]; IC = [5,.1]; R = 1;

x1 = chebfun(5,dom);
x2 = chebfun(.1,dom);
x3 = chebfun(1,dom);
init = [x1,x2];

g = @(x) [IC(1) + cumsum(-R*x(:,1).*x(:,2)), ...
    IC(2) + cumsum(R*x(:,1).*x(:,2)-x(:,2))];

soln = AA(g,init);

%%
F = []; X = []; errVec = []; dom = x.domain;
tol = 1e-8; err = 2*tol; k = 1;

x = stk(x); f = @(x) stk(g(unstk(x,dom))) - x;
x(:,2) = stk(g(unstk(x(:,1),dom)));
errVec(k) = norm(x(:,1)-x(:,2)); k = k+1;
while err > tol
    fk = f(x(:,k));
    X = [ X, x(:,k)-x(:,k-1) ];
    F = [ F, fk - f(x(:,k-1)) ];
    gamma = F \ fk;
    x(:,k+1) = x(:,k) + fk - (X+F)*gamma;

    err = norm(x(:,k)-x(:,k+1)); errVec(k) = err; k = k+1;    
end
soln = unstk(x(:,end),dom); 
numIter = k-1;


%% State Dependance
hist1 = chebfun(@(t)log(t),[0,.1]);
hist2 = chebfun(@(t)1./t,[0,.1]);
g = [join
init1 = chebfun(@(t)log(.1),[.1,1]);
init2 = chebfun(@(t)1./.1,[.1,1]);

g = join