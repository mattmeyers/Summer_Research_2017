function [soln,x,numIter,errVec] = AA(g,x,sol,tol)
%AA     Performs Anderson Acceleration on a scalar ODE or DDE
% Input:
%   g       chebfun being accelerated
%   x       initial guess
%   mk      number of previous terms used, ignore to use full history
% Output:
%   soln    the chebfun solution, last column of x
%   x       all iterations of the solution
%   numIter number of iterations 
%   errVec  vector of successive errors between steps   

% if ~exist('mk','var') || isempty(arg1)
%    mk = inf; 
% end
if ~exist('tol','var')
    tol = 1e-8;
end

mk = inf; %t = sol(:,1); y = sol(:,2);
F = []; X = []; errVec = []; dom = x.domain; 
err = 2*tol; k = 1;

if size(x,2) == 1    % If scalar
    f = @(x) g(x) - x;
    x(:,2) = g(x(:,1)); k = k+1; 
    %errVec = norm(x(:,2)-x(:,1));
    errVec = norm(x(:,2)-sol);
    while err > tol
        if k == 2
            m = min(k-1,mk);
            fk = f(x(:,k));
            X = [ x(:,k)-x(:,k-1) ];
            F = [ fk - f(x(:,k-1)) ];        
        else
            m = min(k-1,mk);
            fk = f(x(:,k));
            X = [ X, x(:,k)-x(:,k-1) ];
            F = [ F, fk - f(x(:,k-1)) ];       
        end
        gamma = F(:,k-m:k-1) \ fk;
        x(:,k+1) = x(:,k) + fk - (X(:,k-m:k-1)+F(:,k-m:k-1))*gamma;
        %err = norm(x(:,k)-x(:,k+1)); errVec(k) = err; k = k+1;
        err = norm(x(:,k+1)-sol); errVec(k) = err; k = k + 1;
    end

    soln = x(:,end);
    numIter = k-1;
else    % If system
    sol = stk(sol);
    x = stk(x); f = @(x) stk(g(unstk(x,dom))) - x;
    x(:,2) = stk(g(unstk(x(:,1),dom)));
    %errVec(k) = norm(x(:,1)-x(:,2)); k = k+1;
    errVec = norm(sol-x(:,2)); k = k+1;
    
    while err > tol
        fk = f(x(:,k));
        X = [ X, x(:,k)-x(:,k-1) ];
        F = [ F, fk - f(x(:,k-1)) ];
        gamma = F \ fk;
        x(:,k+1) = x(:,k) + fk - (X+F)*gamma;
    
        %err = norm(x(:,k)-x(:,k+1)); errVec(k) = err; k = k+1;   
        err = norm(sol-x(:,k+1)); errVec(k) = err; k = k+1;
    end
    soln = unstk(x(:,end),dom); 
    numIter = k-1;
end

end