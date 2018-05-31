%% Solving y'(t)=-0.2*y(t-1)
clear all, close all
yhist = chebfun(@(t)(-1).^floor(-4*t),[-1,0],'splitting','on');
yinit = chebfun(@(t) 1,[0,4]);
y = join(yhist,yinit);

tol = 1e-8; err = 1; numIter = 0; errArr = [];

while err > tol
    y0 = yhist(0);
    yNew = y0 - (1/5)*cumsum(chebfun(@(s)y(s-1),[0,4],'splitting','on'));
    yNew = join(yhist,yNew);
    err = norm(yNew-y);
    errArr = [errArr; err];
    y = yNew;
    numIter = numIter+1;
end

plot(y)
title(sprintf('y''(t)=-0.2*y(t-1): Converged in %d iterations',numIter))
xlabel('t'), ylabel('y(t)')
ylim([-1.5 1.5])
%print -depsc Eq1Soln.eps

%%
% Looking at the error between iterations
semilogy(errArr,'o-')
title('Error Between Iterations'), xlabel('Iteration'), ylabel('Error')
%print -depsc Eq1Err.eps

%% Solving y'(x)=-y(x-1)
clear all, close all
yhist = chebfun(@(t) 1,[-1 0]);
yinit = chebfun(@(t) 1,[0 5]);
y = join(yhist,yinit);

tol = 1e-14; err = 1; numIter = 0; errArr = [];

while err > tol
   y0 = yhist(0);
   yNew = y0 - cumsum(chebfun(@(s)y(s-1),[0,5],'splitting','on'));
   yNew = join(yhist,yNew);
   err = norm(yNew-y);
   errArr = [errArr; err];
   y = yNew;
   numIter = numIter+1;
end

plot(y)
title(sprintf('y''(t)=-y(t-1): Converged in %d iterations',numIter))
xlabel('t'), ylabel('y(t)')
ylim([-.75 1.25])
%print -depsc Eq2Soln.eps

%%
% Looking at the error
semilogy(errArr,'o-')
title('Error Between Iterations'), xlabel('Iteration'), ylabel('Error')
%print -depsc Eq2Err.eps

%% Solving y'(t)=-.5y(t-1)+.2y(t-.25)
clear all, close all
yhist = chebfun(@(t) t.^2,[-1 0]);
yinit = chebfun(@(t) 1,[0 5]);
y = join(yhist,yinit);

tol = 1e-14; err = 1; numIter = 0; errArr = [];

while err > tol
   y0 = yhist(0);
   yNew = y0 + cumsum(chebfun(@(s)-.5*y(s-1)+.2*y(s-.25),[0,5],'splitting','on'));
   yNew = join(yhist,yNew);
   err = norm(yNew-y);
   errArr = [errArr; err];
   y = yNew;
   numIter = numIter+1;
end

plot(y)
figure
semilogy(errArr,'o-')

%% Solving y'(x) = -y(x-1) with ddecd
clear all, close all
ddefun = @(x,y) -.2*y(x-1);
dom = [0,4];
hist = chebfun(@(x)(-1).^floor(-4*x),[-1,0],'splitting','on');
init = chebfun(@(x) 1, dom,'splitting','on');

y = ddecd(ddefun,dom,hist,init);
plot(y)


