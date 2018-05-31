%% A new state dependent DDE, Karoui (2.11), PI
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

%% Comparing Errors
semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o')

%% ddesd
hist = @(t) 1;
delay = @(t,y) t - log(t) - 1;
fun = @(t,y,Z) ((t-1)./t).*y.*Z;

sol = ddesd(fun,delay,hist,dom);


%% Karoui (2.12), PI
clear all, close all
dom = [1,7]; hdom = [0,1]; y0 = 1; 
tol = 1e-8; err = 2*tol; k = 1; errVec = [];
opt = {'splitting','on'};
hist = chebfun(1,hdom,opt{:});
init = chebfun(1,dom,opt{:});
y = join(hist,init);

z = y0 + cumsum(chebfun(@(s) ((s-1)./s).*y(s).*y(s-log(s)-1),dom,opt{:}));
y(:,2) = join(hist,z);
err = norm(y(:,2)-y(:,1)); errVec = err; k = k+1;

while err > tol
    z = y0 + cumsum(chebfun(@(s) ((s-1)./s).*y(s,k).*y(s-log(s)-1,k),dom,opt{:}));
    y(:,k+1) = join(hist,z);
    err = norm(y(:,k+1)-y(:,k)); errVec(k) = err; k = k+1;
end


%% Bellen (16), PI
clear all, close all
dom = [0,4]; hdom = [-1,0]; y0 = 1; 
tol = 1e-8; err = 2*tol; k = 1; errVec = [];
opt = {'splitting','on'};
hist = chebfun(@(s)1-s,hdom,opt{:});
init = chebfun(1,dom,opt{:});
y = join(hist,init);

z = y0 + cumsum(chebfun(@(s) ((s-1)./s).*y(s).*y(s-log(s)-1),dom,opt{:}));
y(:,2) = join(hist,z);
err = norm(y(:,2)-y(:,1)); errVec = err; k = k+1;

while err > tol
    z = y0 + cumsum(chebfun(@(s) ((s-1)./s).*y(s,k).*y(s-log(s)-1,k),dom,opt{:}));
    y(:,k+1) = join(hist,z);
    err = norm(y(:,k+1)-y(:,k)); errVec(k) = err; k = k+1;
end