%% Eq1 y'+abs(y)y=0
clear all, close all
tols = [1e-5,1e-6,1e-7,1e-8,1e-9,1e-10];
%% true solution
dom = [0,5]; y0 = 1;
g = @(y) y0 + cumsum(-abs(y).*y);
init = chebfun(@(t)exp(-t),dom);
[tSoln] = AA(g,init);

%% PI
for i = 1:length(tols)
    tic
    tol = tols(i); err = 2*tol; k1 = 1; errVec1 = [];
    y = chebfun(@(t)1,dom,'splitting','on');
    
    y(:,2) = y0 - cumsum(abs(y).*y);
    err = norm(tSoln-y(:,2)); errVec1 = err; k1 = k1+1;
    
    while err > tol
        y(:,k1+1) = y0 - cumsum(abs(y(:,k1)).*y(:,k1));
        err = norm(tSoln-y(:,k1+1)); errVec1(k1) = err; k1 = k1+1;
    end
    soln1 = y(:,end);
    t1(i) = toc;
end

%% PIAA 
for i = 1:length(tols)
    tic
    g = @(y) y0 + cumsum(-abs(y).*y);
    init = chebfun(@(t)exp(-t),dom);
    [soln2,~,k2,errVec2] = AA(g,init,tSoln,tols(i));
    t2(i) = toc;
end


%% Newton
for i = 1:1;length(tols)
    tic
    dom = [0,5]; tol = tols(i); err = 2*tol; k3 = 1; errVec3 = [];
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
    t3(i) = toc;
end

%%
clf
loglog(t1,tols,'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','g'), hold on
loglog(t2,tols,'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','b')
loglog(t3,tols(1:5),'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','r')
xlabel('time (s)'), ylabel('tolerance'), legend('PI','PIAA','Newton')
fig = gcf; fig.PaperUnits='inches'; fig.PaperPosition=[0,0,7.5,4.43];



%% Eq3: y''+sin(y)=0
clear all, close all
%% Chebfun solution
N = chebop(@(t,y)diff(y,2)+sin(y),[0,5]);
N.lbc = [1;0];
tSoln = N\0;
tSoln = [tSoln,diff(tSoln)];

%% PI

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
