%% Applying Newton to y'(t)=-y(t)y(t-1)
clear all, close all
lag = @(f,d) newDomain(f,[f.domain(1)+d,f.domain(end)+d]);
dom = [0,2]; hdom = [-1,0]; opt = {'splitting','on'};
hist  = chebfun(1,hdom,'splitting','on'); y0 = 1;
y = join(hist,chebfun(1,dom,'splitting','on'));

y1 = restrict(y,dom);
y2 = restrict(lag(y,1),dom);
%restrict(lag(join(chebfun(0,hdom),d),1),dom)

N = chebop(@(t,d) d + ddeVolt(@(t,s)y1(s),d) + volt(@(t,s)y2(s),d),dom);
Z = y0 - restrict(y,dom) - cumsum(chebfun(@(s)y(s-1).*y(s),dom,'splitting','on'));

delta = N\Z;
yn = y + join(chebfun(0,hdom),delta);

%%
tol = 1e-8; err = 2*tol; k = 1; errVec1 = [];
dom = [0,2]; hdom = [-1,0]; y0 = -1;
yhist = chebfun(-1,hdom); yinit = chebfun(-1,dom);
y = join(yhist,yinit);

z = y0 - cumsum(chebfun(@(s)y(s).*y(s-1),dom,'splitting','on')); 
y(:,2) = join(yhist,z);
err = norm(y(:,1)-y(:,2)); errVec1 = err; k = k+1;

while err > tol
   z = y0 - cumsum(chebfun(@(s)y(s,k).*y(s-1,k),dom,'splitting','on'));
   y(:,k+1) = join(yhist,z);
   err = norm(y(:,k)-y(:,k+1)); errVec1(k) = err; k = k+1;
end
plot(y(:,end))


%% State dependent system PIAA
clear all, close all
hdom = [1e-5,.1]; dom = [.1,1]; y10 = log(.1); y20 = 1/.1;
hist1 = chebfun(@(t)log(t),hdom); hist2 = chebfun(@(t)1./t,hdom);
g = @(y) [join(hist1,y10+cumsum(y(:,1))),join(hist2,y20-cumsum(y(exp(1-y(:,2)),2).*y(:,2).^2.*exp(1-y(:,2))))];
init = [chebfun(y10,dom),chebfun(y20,dom)];
[soln,~,numIter,errVec] = AA(g,init);