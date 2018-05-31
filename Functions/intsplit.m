function soln = intsplit(g,x,ic)
%INTSPLIT   Solves an ODE by repeated integration
% Inputs:
%   g       Function, y^(n)(x)=g(x,y,y',...,y^(n-1))
%   x       Initial chebfun guess
%   ic      Initial conditions, ordered in decreasing derivative order
% Output:
%   soln    Final solution, chebfun

soln = intsplit_rec(g,x,ic);
end

function [soln,ic] = intsplit_rec(g,x,ic)
[a,b] = domain(x);
[y,ic] = splitsolve(g,x,ic);
if ~(isa(y,'chebfun'))
    [x1,ic] = intsplit_rec(g,restrict(x,[a,a+(b-a)/2]),ic);
    [x2,ic] = intsplit_rec(g,restrict(x,[a+(b-a)/2,b]),ic);
    soln = join(x1,x2);
else
    soln = y;
end

end

function [soln,ic] = splitsolve(g,x,ic)
tol = 1e-8; err = 2*tol; k = 1;
while err > tol
   z = ic(1) + cumsum(g(x(:,k)));
   for i = 2:length(ic)
       z(:,i) = ic(i) + cumsum(z(:,i-1));
   end
   x(:,k+1) = z(:,end);
   
   if norm(x(:,k+1)) > 1e5
      soln = 0;
      return; 
   end
   err = norm(x(:,k+1)-x(:,k)); k = k+1; 
end
soln = x(:,end);
for i = 1:length(ic)
    ic(i) = z(end,i);
end

end

