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