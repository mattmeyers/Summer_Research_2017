%% Solving the Kermack-McKendrick Equations
close all, clear all
hdom = [-1 0]; dom = [0 4];

hist1 = chebfun(@(t) 5, hdom, 'splitting','on');
hist2 = chebfun(@(t) .1, hdom, 'splitting','on');
hist3 = chebfun(@(t) 1, hdom, 'splitting','on');
hist = [hist1 hist2 hist3];

init1 = chebfun(@(t) 5, dom, 'splitting','on');
init2 = chebfun(@(t) .1, dom, 'splitting','on');
init3 = chebfun(@(t) 1, dom, 'splitting','on');
init = [init1 init2 init3];

y = [join(hist1,init1) join(hist2,init2) join(hist3,init3)];
y0 = [init1(0) init2(0) init3(0)];

err = 1; tol = 1e-8; numIter = 0;

while err > tol
    y1 = y(:,1); y2 = y(:,2); y3 = y(:,3);
    y1N = y0(1) + cumsum(chebfun(@(s) -y1(s)*y2(s-.1)+y2(s-1),dom,'splitting','on'));
    y2N = y0(2) + cumsum(chebfun(@(s) y1(s)*y2(s-.1)-y2(s),dom,'splitting','on'));
    y3N = y0(3) + cumsum(chebfun(@(s) y2(s)-y2(s-1),dom,'splitting','on'));
    
    y1N = join(hist1,y1N); y2N = join(hist2,y2N); y3N = join(hist3,y3N);
    err = [norm(y1-y1N) norm(y2-y2N) norm(y3-y3N)]; err = max(err);
    
    y = [y1N y2N y3N]; numIter = numIter + 1;
end

plot(y(:,1)), hold on, plot(y(:,2)), plot(y(:,3))

%% Solving a different type of DDE
clear all, close all
dom = [0 50]; err = 1; tol = 1e-12; numIter = 0; errArr = [];

y0 = 0.5; y = chebfun(@(s) .5,dom,'splitting','on');

while err > tol
   yN = y0 + cumsum(chebfun(@(s) -.25*y(s)+y(s/2)*(1-y(s/2)),dom));
   err = norm(yN-y);
   errArr = [errArr err];
   y = yN;
   numIter = numIter + 1;
end

plot(y)
title('Eq 2 Solution'), xlabel('t'), ylabel('y(t)')
%print -depsc Eq2Soln.eps
figure
semilogy(errArr,'o-')
title('Eq 2 Error'), xlabel('Iteration'), ylabel('Error')
%print -depsc Eq2Err.eps



