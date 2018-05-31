%% Solving Kermack McKendrick ODE without acceleration
clear all, close all
dom = [0,10];

x1 = chebfun(5,dom);
x2 = chebfun(.1,dom);
x3 = chebfun(1,dom);

g1 = @(x1,x2,x3) x1(0) + cumsum(-x1.*x2);
g2 = @(x1,x2,x3) x2(0) + cumsum(x1.*x2-x2);
g3 = @(x1,x2,x3) x3(0) + cumsum(x2);

for k = 1:10
    x1(:,k+1) = g1(x1(:,k),x2(:,k),x3(:,k));
    x2(:,k+1) = g2(x1(:,k),x2(:,k),x3(:,k));
    x3(:,k+1) = g3(x1(:,k),x2(:,k),x3(:,k));
end

%% Solving Kermack McKendrick with acceleration
clear all, close all
dom = [0,18];

x1 = chebfun(@(x) 1./x,dom);
x2 = chebfun(@(x)exp(-(x-7).^2),dom);
x3 = chebfun(@(x)log(x),dom);

g1 = @(x1,x2,x3) x1(0) + cumsum(-2*x1.*x2);
g2 = @(x1,x2,x3) x2(0) + cumsum(2*x1.*x2-x2);
g3 = @(x1,x2,x3) x3(0) + cumsum(x2);

f1 = @(x1,x2,x3) g1(x1,x2,x3) - x1;
f2 = @(x1,x2,x3) g2(x1,x2,x3) - x2;
f3 = @(x1,x2,x3) g3(x1,x2,x3) - x3;
F1 = []; F2 = []; F3 = [];
X1 = []; X2 = []; X3 = [];
x1(:,2) = g1(x1(:,1),x2(:,1),x3(:,1));
x2(:,2) = g2(x1(:,1),x2(:,1),x3(:,1));
x3(:,2) = g3(x1(:,1),x2(:,1),x3(:,1));

for k = 2:10
    m = min(k-1,6);
    X1 = [ X1, x1(:,k)-x1(:,k-1) ];
    X2 = [ X2, x2(:,k)-x2(:,k-1) ];
    X3 = [ X3, x3(:,k)-x3(:,k-1) ];
    
    F1 = [ F1, f1(x1(:,k),x2(:,k),x3(:,k)) - f1(x1(:,k-1),x2(:,k-1),x3(:,k-1)) ];
    F2 = [ F2, f2(x1(:,k),x2(:,k),x3(:,k)) - f2(x1(:,k-1),x2(:,k-1),x3(:,k-1)) ];
    F3 = [ F3, f3(x1(:,k),x2(:,k),x3(:,k)) - f3(x1(:,k-1),x2(:,k-1),x3(:,k-1)) ];
    
    gamma1 = F1 \ f1(x1(:,k),x2(:,k),x3(:,k));
    gamma2 = F2 \ f2(x1(:,k),x2(:,k),x3(:,k));
    gamma3 = F3 \ f3(x1(:,k),x2(:,k),x3(:,k));
    
    x1(:,k+1) = x1(:,k) + f1(x1(:,k),x2(:,k),x3(:,k)) - (X1+F1)*gamma1;
    x2(:,k+1) = x2(:,k) + f2(x1(:,k),x2(:,k),x3(:,k)) - (X2+F2)*gamma2;
    x3(:,k+1) = x3(:,k) + f3(x1(:,k),x2(:,k),x3(:,k)) - (X3+F3)*gamma3;
end

plot(x1(:,end)), hold on, plot(x2(:,end)), plot(x3(:,end))