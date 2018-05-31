%% Eq1 y'+abs(y)y=0
clear all, close all

%% true solution
dom = [0,5]; y0 = 1;
g = @(y) y0 + cumsum(-abs(y).*y);
init = chebfun(@(t)exp(-t),dom);
[tSoln] = AA(g,init);

%% PI
tic
tol = 1e-8; err = 2*tol; k1 = 1; errVec1 = [];
y = chebfun(@(t)1,dom,'splitting','on');

y(:,2) = y0 - cumsum(abs(y).*y); 
err = norm(tSoln-y(:,2)); errVec1 = err; k1 = k1+1;

while err > tol
   y(:,k1+1) = y0 - cumsum(abs(y(:,k1)).*y(:,k1)); 
   err = norm(tSoln-y(:,k1+1)); errVec1(k1) = err; k1 = k1+1;
end
soln1 = y(:,end);
t1 = toc;

%% PIAA 
tols = [1e-5,1e-6,1e-7,1e-8,1e-9,1e-10]; times = [];
for i = 1:length(tols)
    tic
    g = @(y) y0 + cumsum(-abs(y).*y);
    init = chebfun(@(t)exp(-t),dom);
    [soln2,~,k2,errVec2] = AA(g,init,tSoln,tols(i));
    times(i) = toc;
end


%% Newton
tic
dom = [0,5]; err = 2*tol; k3 = 1; errVec3 = [];
y0 = 1; yh = chebfun(@(t)1,dom);
N = chebop(@(t,d)d+volt(@(t,s)abs(yh(t)),d)+volt(@(t,s)yh(t).*sign(yh(t)),d),dom);
r = y0-cumsum(abs(yh).*yh)-yh;

y = N\r;
y = yh+y;
errVec3 = norm(tSoln-y); k3 = k3 + 1;

while err > tol
    yh = y(:,1);
    N = chebop(@(t,d)d+volt(@(t,s)abs(yh(t)),d)+volt(@(t,s)yh(t).*sign(yh(t)),d),dom);
    r = y0-cumsum(abs(yh).*yh)-yh;

    yn = N\r;
    yn = yh+yn;
    err = norm(tSoln-yn); errVec3(k3+1) = err; y = yn; k3 = k3+1;
end
soln3 = y;
t3 = toc;

%% Comparing Error
clf
semilogy(errVec1,'-s','linewidth',2,'markersize',7), hold on
semilogy(errVec2,'-s','linewidth',2,'markersize',7)
semilogy(errVec3,'-s','linewidth',2,'markersize',7)
xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA','Newton'), ylim([1e-10,1e4])


%% Eq2 y' + y^2 = 0
clear all, close all
%% True solution
dom = [0,2]; y0 = 1;
N = chebop(@(t,y)diff(y)+y.^2,dom);
N.lbc = y0;
tSoln = N\0;

%% PI
tic
tol = 1e-8; err = 2*tol; k1 = 1; errVec1 = [];
y = chebfun(@(t)1,dom,'splitting','on');

y(:,2) = y0 - cumsum(y.^2); 
err = norm(y(:,1)-y(:,2)); errVec1 = err; k1 = k1+1;

while err > tol
   y(:,k1+1) = y0 - cumsum(y(:,k1).^2); 
   err = norm(y(:,k1)-y(:,k1+1)); errVec1(k1) = err; k1 = k1+1;
end
soln1 = y(:,end);
t1 = toc;

%% PIAA 
tic
g = @(y) y0 + cumsum(-y.^2);
init = chebfun(@(t)exp(-t),dom);
[soln2,~,k2,errVec2] = AA(g,init,tSoln);
t2 = toc;

%% Newton
tic
err = 2*tol; k3 = 1; errVec3 = [];
y0 = 1; yh = chebfun(@(t)1,dom);
N = chebop(@(t,d)d+2*volt(@(t,s)yh(t),d),dom);
r = y0-cumsum(abs(yh).*yh)-yh;

y = N\r;
y = yh+y;
errVec3 = norm(tSoln-y); k3 = k3 + 1;

while err > tol
    yh = y(:,1);
    N = chebop(@(t,d)d+2*volt(@(t,s)yh(t),d),dom);
    r = y0-cumsum(abs(yh).*yh)-yh;

    yn = N\r;
    yn = yh+yn;
    err = norm(tSoln-yn); errVec3(k3) = err; y = yn; k3 = k3+1;
end
soln3 = y;
t3 = toc;

%% Comparing Error
semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o'),semilogy(errVec3,'-o')
xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA','Newton')



%% Eq3: y''+sin(y)=0
clear all, close all
%% Chebfun solution
N = chebop(@(t,y)diff(y,2)+sin(y),[0,5]);
N.lbc = [1;0];
tSoln = N\0;
tSoln = [tSoln,diff(tSoln)];

%% PI
dom = [0,5];
u = chebfun(@(t) 1, dom); v = chebfun(@(t) 0, dom); y = [u,v];
tol = 1e-8; err = 2*tol; k1 = 1; errVec1 = [];

while err > tol
    u = 1 + cumsum(chebfun(@(t) y(t,2),dom));
    v = cumsum(chebfun(@(t) -sin(y(t,1)),dom));
    yn = [u,v];
    err = norm(yn-y); errVec1(k1) = err; y = yn;
    k1 = k1+1;
end
soln1 = y;

%% PIAA
g = @(y)[1+cumsum(y(:,2)),cumsum(-sin(y(:,1)))];
init = [chebfun(1,dom),chebfun(0,dom)];
[soln2,~,k2,errVec2] = AA(g,init,tSoln);

%% Newton
dom = [0,5]; tol = 1e-8; err = 2*tol; k3 = 1; errVec3 = [];
y0 = 1; y1h = chebfun(@(t)1,dom); y2h = chebfun(@(t)0,dom); yh = [y1h,y2h];
N = chebop(@(t,d,e)[d-volt(@(t,s)1,e);e+volt(@(t,s)cos(y1h(t)),d)],dom);
r = [1+cumsum(y2h)-y1h;-cumsum(sin(y1h))-y2h];

y = N\r;
y = [y1h+y{1},y2h+y{2}];
errVec3 = norm(tSoln-y); k3 = k3 + 1;

while err > tol
    y1h = y(:,1); y2h = y(:,2); yh = [y1h,y2h];
    N = chebop(@(t,d,e)[d-volt(@(t,s)1,e);e+volt(@(t,s)cos(y1h(t)),d)],dom);
    r = [1+cumsum(y2h)-y1h;-cumsum(sin(y1h))-y2h];
    
    yn = N\r;
    yn = [y1h+yn{1},y2h+yn{2}];
    err = norm(tSoln-yn); errVec3(k3) = err; y = yn; k3 = k3+1;
end
soln3 = y;

%% Comparing Error
clf
semilogy(errVec1,'-s','linewidth',2,'markersize',7), hold on
semilogy(errVec2,'-s','linewidth',2,'markersize',7)
semilogy(errVec3,'-s','linewidth',2,'markersize',7)
xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA','Newton')


%% Eq4: y'+y^3=3cos(t)
clear all, close all
dom = [0,20]; y0 = 0; tol = 1e-8;
%% True solution
N = chebop(@(t,y) diff(y)+y^3,dom);
N.lbc = y0;
t = chebfun('t',dom);
tSoln = N\(3*cos(t));

%% PI
err = 2*tol; k1 = 1; errVec1 = [];
y = chebfun(@(t)1,dom); t = chebfun('t',dom);

y(:,2) = y0 - cumsum(y.^3)+cumsum(3*cos(t)); 
err = norm(tSoln-y(:,2)); errVec1 = err; k1 = k1+1;

while err > tol
   y(:,k1+1) = y0 - cumsum(y(:,k1).^3)+cumsum(3*cos(t)); 
   err = norm(tSoln-y(:,k1+1)); errVec1(k1) = err; k1 = k1+1;
end
soln1 = y(:,end);

%% PIAA
g = @(y) y0 - cumsum(y.^3)+cumsum(3*cos(t));
init = chebfun(@(t)y0,dom);
[soln2,~,k2,errVec2] = AA(g,init,tSoln);

%% Newton
err = 2*tol; k3 = 1; errVec3 = [];
yh = chebfun(@(t)0,dom);
N = chebop(@(t,d)d+3*volt(@(t,s)yh(s).^2,d),dom);
r = y0-cumsum(yh.^3)+cumsum(chebfun(@(t)3*cos(t),dom))-yh;

y = N\r;
y = yh+y;
errVec3 = norm(tSoln-y); k3 = k3 + 1;

while err > tol
    yh = y(:,1);
    N = chebop(@(t,d)d+3*volt(@(t,s)yh(s).^2,d),dom);
    r = y0-cumsum(yh.^3)+cumsum(chebfun(@(t)3*cos(t),dom))-yh;

    yn = N\r;
    yn = yh+yn;
    err = norm(tSoln-yn); errVec3(k3) = err; y = yn; k3 = k3+1;
end
soln3 = y;

%% Comparing Error
clf
semilogy(errVec1,'-s','linewidth',2,'markersize',7), hold on
semilogy(errVec2,'-s','linewidth',2,'markersize',7)
semilogy(errVec3,'-s','linewidth',2,'markersize',7)
xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA','Newton')



%% Eq5: y' = -.2*y(t-1), phi = (-1)^floor(-4t)
clear all, close all
hdom = [-1,0]; dom = [0,4]; y0 = 1; tol = 1e-8;

%% True solution
hist = @(t) (-1)^floor(-4*t);
delay = @(t,y) t-1;
fun = @(t,y,Z) -.2*Z;
sol = ddesd(fun,delay,hist,dom);
g = @(t)deval(sol,t);
tSoln = chebfun(g,dom,'splitting','on');

%% true solution 2
yhist = chebfun(@(t)(-1).^floor(-4*t),hdom,'splitting','on');
yinit = chebfun(@(t) 1,dom);
y = join(yhist,yinit);

err = 2*tol; k1 = 1; errVec1 = [];

while err > 1e-12
    z = y0 - (1/5)*cumsum(chebfun(@(s)y(s-1,k1),dom,'splitting','on'));
    y(:,k1+1) = join(yhist,z);
    err = norm(y(:,k1)-y(:,k1+1));
    errVec1(k1) = err; k1 = k1+1;
end
tSoln = y(:,end);

%% PI
yhist = chebfun(@(t)(-1).^floor(-4*t),hdom,'splitting','on');
yinit = chebfun(@(t) 1,dom);
y = join(yhist,yinit);

err = 2*tol; k1 = 1; errVec1 = [];

while err > tol
    z = y0 - (1/5)*cumsum(chebfun(@(s)y(s-1,k1),dom,'splitting','on'));
    y(:,k1+1) = join(yhist,z);
    err = norm(tSoln-y(:,k1+1));
    errVec1(k1) = err; k1 = k1+1;
end
soln1 = y(:,end);

%% PIAA
g = @(y) join(yhist,y0-.2*cumsum(chebfun(@(s)y(s-1),dom,'splitting','on')));
init = join(yhist,yinit);
[soln2,~,k2,errVec2] = AA(g,init,tSoln);

%% Comparing error
clf
semilogy(errVec1,'-s','linewidth',2,'markersize',7), hold on
semilogy(errVec2,'-s','linewidth',2,'markersize',7)
xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA')



%% Eq6: y' = -.25*y+y(t/2)(1-y(t/2))
clear all, close all
dom = [0,25]; y0 = .5; tol = 1e-8;

%% True Solution
dom = [0 25]; err = 1; tol = 1e-12; numIter = 0; errArr = [];

y0 = 0.5; y = chebfun(@(s) .5,dom,'splitting','on');

while err > tol
   yN = y0 + cumsum(chebfun(@(s) -.25*y(s)+y(s/2)*(1-y(s/2)),dom));
   err = norm(yN-y);
   errArr = [errArr err];
   y = yN;
   numIter = numIter + 1;
end
tSoln = y(:,end);

%% PI
dom = [0 25]; err = 1; tol = 1e-8; k1 = 1; errVec1 = [];

y0 = 0.5; y = chebfun(@(s) .5,dom,'splitting','on');

while err > tol
   yN = y0 + cumsum(chebfun(@(s) -.25*y(s)+y(s/2)*(1-y(s/2)),dom));
   err = norm(yN-tSoln);
   errVec1(k1) = err;
   y = yN;
   k1 = k1 + 1;
end
soln1 = y(:,end);

%% PIAA
g = @(y) y0+cumsum(chebfun(@(s)-.25*y(s)+y(s/2).*(1-y(s/2)),dom,'splitting','on'));
init = chebfun(.5,dom,'splitting','on');
[soln2,~,k2,errVec2] = AA(g,init,tSoln);

%% Comparing error
clf
semilogy(errVec1,'-s','linewidth',2,'markersize',7), hold on
semilogy(errVec2,'-s','linewidth',2,'markersize',7)
xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA')



%% Eq7: y' = (1/t)*y*y(ln(y))
clear all, close all
%% True solution
dom = [1,5]; hdom = [0,1]; y0 = 1; 
tol = 1e-12; err = 2*tol; k = 1;
opt = {'splitting','on'};
hist = chebfun(1,hdom,opt{:});
init = chebfun(1,dom,opt{:});
y = join(hist,init);

z = y0 + cumsum(chebfun(@(s) (1./s).*y(s).*y(log(y(s))),dom,opt{:}));
y(:,2) = join(hist,z);
err = norm(y(:,2)-y(:,1)); k = k+1;

while err > tol
    z = y0 + cumsum(chebfun(@(s) (1./s).*y(s,k).*y(log(y(s,k)),k),dom,opt{:}));
    y(:,k+1) = join(hist,z);
    err = norm(y(:,k+1)-y(:,k)); k = k+1;
end
tSoln = y(:,end);

%% PI
dom = [1,5]; hdom = [0,1]; y0 = 1; 
tol = 1e-8; err = 2*tol; k1 = 1; errVec1 = [];
opt = {'splitting','on'};
hist = chebfun(1,hdom,opt{:});
init = chebfun(1,dom,opt{:});
y = join(hist,init);

z = y0 + cumsum(chebfun(@(s) (1./s).*y(s).*y(log(y(s))),dom,opt{:}));
y(:,2) = join(hist,z);
err = norm(y(:,2)-tSoln); errVec1 = err; k1 = k1+1;

while err > tol
    z = y0 + cumsum(chebfun(@(s) (1./s).*y(s,k1).*y(log(y(s,k1)),k1),dom,opt{:}));
    y(:,k1+1) = join(hist,z);
    err = norm(y(:,k1+1)-tSoln); errVec1(k1) = err; k1 = k1+1;
end
soln1 = y(:,end);

%% PIAA
g = @(y) join(hist, y0 + cumsum(chebfun(@(s) (1./s).*y(s).*y(log(y(s))),dom,'splitting','on')));
in = join(hist,init);
[soln2,~,k2,errVec2] = AA(g,in,tSoln);

%% Comparing error
clf
semilogy(errVec1,'-s','linewidth',2,'markersize',7), hold on
semilogy(errVec2,'-s','linewidth',2,'markersize',7)
xlabel('Iteration'),ylabel('Error')
legend('PI','PIAA')