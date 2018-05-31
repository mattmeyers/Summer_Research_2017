x = 1.29;
g = @(x) -x^2 + 5*x - 3.5;
for k = 1:20
    x(k+1) = g(x(k));
end

%%
f = @(x) g(x) - x;
F = [];
X = [];
x(2) = g(x(1));
for k = 2:10
    m = k;
    X = [ X, x(k)-x(k-1) ];
    F = [ F, f(x(k)) - f(x(k-1)) ];
    gamma = F \ f(x(k));
    x(k+1) = x(k) + f(x(k)) - (X+F)*gamma;
end

%%
x = chebfun(1,[0 5]);
g = @(x) x(0) + cumsum(exp(x/10));
for k = 1:10
    x(:,k+1) = g(x(:,k));
end

%%
f = @(x) g(x) - x;
F = [];
X = [];
x(:,2) = g(x(:,1));
for k = 2:10
    m = min(k-1,6);
    fk = f(x(:,k));
    X = [ X(:,1:m-1), x(:,k)-x(:,k-1) ];
    F = [ F(:,1:m-1), fk - f(x(:,k-1)) ];
    gamma = F \ fk;
    x(:,k+1) = x(:,k) + fk - (X+F)*gamma;
end