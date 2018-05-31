%% Solving a system of ODEs without Chebfun
L = chebop(0, 5);
L.op = @(x, u, v) [diff(u) - v; diff(v) + u];
L.lbc = @(u, v) [u-1; v];
rhs = [0; 0];
U = L\rhs;
plot(U)

%% Solving with Picard, no AA
clear all, close all
dom = [0,5];
u = chebfun(@(t) 1, dom); v = chebfun(@(t) 0, dom);
tol = 1e-8; err = 2*tol; k = 1; errVec1 = [];

while err > tol
    u(:,k+1) = 1 + cumsum(chebfun(@(t) v(t,k),dom));
    v(:,k+1) = cumsum(chebfun(@(t) -u(t,k),dom));
    err = max(norm(u(:,k+1)-u(:,k)),norm(v(:,k+1)-v(:,k))); errVec1(k) = err;
    k = k+1;
end

plot(u(:,end)), hold on, plot(v(:,end))

%% Solving with Picard, with AA
dom = [0,5];
g = @(x) [1 + cumsum(chebfun(@(t) x(t,2),dom,'splitting','on')),cumsum(chebfun(@(t) -x(t,1),dom,'splitting','on'))];
x = [chebfun(@(t) 1, dom,'splitting','on'), chebfun(@(t) 0, dom,'splitting','on')];
x = stk(x);
f = @(x) stk(g(unstk(x,dom)))-x;
F = []; X = []; tol = 1e-8; err = 2*tol; errVec2 = []; k = 1;

x(:,2) = stk(g(unstk(x(:,1),dom)));
errVec2(k) = norm(x(:,1)-x(:,2)); k = k+1;
while err > tol
    fk = f(x(:,k));
    X = [ X, x(:,k)-x(:,k-1) ];
    F = [ F, fk - f(x(:,k-1)) ];
    gamma = F \ fk;
    x(:,k+1) = x(:,k) + fk - (X+F)*gamma;
    
    err = norm(x(:,k)-x(:,k+1)); errVec2(k) = err; k = k+1;    
end
soln = unstk(x(:,end),dom);
figure, plot(soln)
figure, semilogy(errVec1,'-o'), hold on, semilogy(errVec2,'-o')