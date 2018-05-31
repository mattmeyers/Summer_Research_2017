%% Solving y'(t)=sin(t^2)y
clear all, close all
y0 = 1; dom = [0,8]; err = 1; tol = 1e-14;

y = chebfun(@(x) 1, dom);
k = 1; errArr = [];

while err > tol
    y(:,k+1) = y0 + cumsum(chebfun(@(t)sin(t^2)*y(t,k),dom));
    err = norm(y(:,k+1)-y(:,k));
    errArr(k) = err;
    k = k+1;
end

plot(y(:,end)), title(sprintf('y''(t)=sin(t^2)y(t): %d iterations',k-1))
xlabel('t'), ylabel('y(t)')
%print -depsc Eq1Soln.eps
figure, semilogy(errArr,'-o')
title(sprintf('Successive Error: tol=%g',tol))
xlabel('iteration'), ylabel('Error')
%print -depsc Eq1Err.eps

%% Repeating with AA
g = @(x) y0 + cumsum(chebfun(@(t)sin(t^2).*x(t),dom));
x = chebfun(@(x) 1,dom);
f = @(x) g(x) - x;
F = []; X = []; errVec = []; err = 1; tol = 1e-8; k = 2;
x(:,2) = g(x(:,1));
while err > tol
    fk = f(x(:,k));
    X = [ X, x(:,k)-x(:,k-1) ];
    F = [ F, fk - f(x(:,k-1)) ];
    gamma = F \ fk;
    x(:,k+1) = x(:,k) + fk - (X+F)*gamma;
    
    err = norm(x(:,k)-x(:,k+1)); errVec(k-1) = err; k = k+1;
end
figure 
plot(x(:,end)), title('y''(t)=sin(t^2)y(t) with AA')
xlabel('t'), ylabel('y(t)')%, print -depsc Eq1AASoln.eps

%% Comparing Errors
figure
semilogy(errArr,'-o'), hold on, semilogy(errVec,'-o')
legend('Without AA','With AA')
title('Comparing Errors'), xlabel('Iteration'), ylabel('Error')
print -depsc Eq1CompErr.eps


%% Solving y'(t)=sign(sin(t^2))y
clear all, close all
y0 = 1; dom = [0,8]; err = 1; tol = 1e-8;

y = chebfun(@(x) 1, dom,'splitting','on');
k = 1; errArr = [];

while err > tol
    y(:,k+1) = y0 + cumsum(chebfun(@(t)sign(sin(t^2))*y(t,k),dom,'splitting','on'));
    err = norm(y(:,k+1)-y(:,k));
    errArr(k) = err;
    k = k+1;
end

plot(y(:,end)), title(sprintf('y''(t)=sign(sin(t^2))y(t): %d iterations',k-1))
xlabel('t'), ylabel('y(t)')
%print -depsc Eq2Soln.eps
figure, semilogy(errArr,'-o')
title(sprintf('Successive Error: tol=%g',tol))
xlabel('iteration'), ylabel('Error')
%print -depsc Eq2Err.eps


%% Solving y'(t)=-|y|^2 * y + 3cos(t)
clear all, close all
y0 = 0; dom = [0,5]; err = 1; tol = 1e-8;

y = chebfun(@(x)sin(x),dom);
k=1; errArr = [];

while err > tol
    y(:,k+1) = y0 + cumsum(chebfun(@(t)-abs(y(t,k))^2*y(t,k)+3*cos(t),dom));
    err = norm(y(:,k+1)-y(:,k));
    errArr(k) = err;
    k = k+1;
end

plot(y(:,end)), title(sprintf('y''(t)=sign(sin(t^2))y(t): %d iterations',k-1))
xlabel('t'), ylabel('y(t)')
%print -depsc Eq3Soln.eps
figure, semilogy(errArr,'-o')
title(sprintf('Successive Error: tol=%g',tol))
xlabel('iteration'), ylabel('Error')
%print -depsc Eq3Err.eps

%% Repeating with AA
clear all, close all
dom = [0,20];
g = @(y) cumsum(chebfun(@(t)-abs(y(t)).^2.*y(t)+3*cos(t),dom,'splitting','on'));
%t = chebfun('t',dom);
%g = @(y) cumsum(-abs(y).^2.*y+3*cos(t));
init = chebfun(@(x)cos(x),dom);

[soln,x,numIter,errVec] = AA(g,init);

plot(soln), title(sprintf('y''(t)=|y(t)|^2y(t)+3cos(t): %d iterations',numIter))
xlabel('t'), ylabel('y(t)')
%print -depsc Eq3Soln.eps
figure, semilogy(errVec,'-o')
title(sprintf('Successive Error: tol=%g',1e-8))
xlabel('iteration'), ylabel('Error')
%print -depsc Eq3Err.eps