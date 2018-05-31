function Script0618fun
% Applying the AA algorithm to expectation maximization algorithm
clear all, close all
alpha = [.3,.3,.4]; sigma = [1,1,1]; N = 10000;
mu1 = [0,2,4]; mu2 = [0,1,2]; mu3 = [0,.5,1];

rng(1234);

mu = mu3;
sig = diag(sigma);
R = chol(sig);
z = repmat(mu,N,1) + randn(N,3)*R;

% EM algorithm without AA
    function sum = p(x)
       sum = zeros(N,1);
       for k = 1:3
          sum = sum + alpha(k)*q(x,k);
       end
    end

    function y = q(x,j)
       y = (1/(sqrt(2*pi)*sigma(j))).*exp((-(x-eta(j)).^2)./(2*sigma(j)^2));
    end

    function next = iter(x,j)
        s = (alpha(j)*q(x,j))./(p(x));
        next = (x'*s)/sum(s);
    end

err = 1; tol = 1e-8; numIter = 0; eta = [1,1,1]; res = [];

while err > tol
    for i = 1:3
        etaN(i) = iter(z(:,i),i);
    end
    err = norm(mu-etaN);
    eta = etaN;
    numIter = numIter + 1;
    res = [res err];
    if numIter > 125
        disp('Breaking')
        break
    end
end
plot(res,'o-')


end