%% Solving the Luggage DDE
close all, clear all
hdom = [-.1 0]; dom = [0 12];
gamma = 2.48; beta = 1; tau = 0.1; A = .75; omega = 1.37;
eta = asin(gamma/A);

hist = chebfun(@(t) 0, hdom, 'splitting','on');

init1 = chebfun(@(t) sin(t), dom, 'splitting','on');
init2 = chebfun(@(t) sin(t), dom, 'splitting','on');
init = [init1 init2];

y = [join(hist,init1) join(hist,init2)];
y0 = [hist(0) hist(0)];

err = 1; tol = 1e-8; numIter = 0;

while err > tol
    y1 = y(:,1); y2 = y(:,2);
    y1N = y0(1) + cumsum(chebfun(@(s) y2(s),dom,'splitting','on'));
    y2N = y0(2) + cumsum(chebfun(@(s) A*sin(omega*s+eta)-sign(y1(s))*gamma*cos(y1(s))+sin(y1(s))-beta*y1(s-tau),dom,'splitting','on'));
    
    y1N = join(hist,y1N); y2N = join(hist,y2N);
    err = [norm(y1-y1N) norm(y2-y2N)]; err = max(err);
    
    y = [y1N y2N]; numIter = numIter + 1;
end

plot(y(:,1)), hold on, plot(y(:,2))
