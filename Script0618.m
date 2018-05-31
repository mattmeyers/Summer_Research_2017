%% Solving the system of state dependent DDEs
clear all, close all
hdom = [1e-4,.1]; dom = [.1,1];

hist1 = chebfun(@(t) log(t),hdom,'splitting','on');
hist2 = chebfun(@(t) 1./t,hdom,'splitting','on');

init1 = chebfun(@(t) log(t+.1),dom,'splitting','on');
init2 = chebfun(@(t) (1./t)+.1,dom,'splitting','on');

y = [join(hist1,init1) join(hist2,init2)];
y0 = [init1(0.1) init2(0.1)];

err = 1; tol = 1e-8; numIter = 0;

while err > tol
    y1N = y0(1) + cumsum(chebfun(@(s) y(s,2),dom,'splitting','on'));
    y2N = y0(2) + cumsum(chebfun(@(s) -y(exp(1-y(s,2)),2).*y(s,2).^2.*exp(1-y(s,2)),dom,'splitting','on'));
    
    y1N = join(hist1,y1N); y2N = join(hist2,y2N);
    err = [norm(y(:,1)-y1N) norm(y(:,2)-y2N)]; err = max(err);
    
    y = [y1N y2N]; numIter = numIter + 1;
end

plot(y(:,1)), hold on, plot(y(:,2))

