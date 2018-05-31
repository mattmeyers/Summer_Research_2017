function [soln,x,k,errVec] = AASys(g,dom,x)
x = stk(x);
f = @(x) stk(g(unstk(x,dom)))-x;
F = []; X = []; tol = 1e-8; err = 2*tol; errVec = []; k = 1; 

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

end