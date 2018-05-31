%% Let's try solving y''+sign(y)=0 again
clear all, close all
dom = [0,20]; tol = 1e-12; err = 2*tol; k = 1;

y = chebfun(@(t)cos(t),dom,'splitting','on');

while err > tol
   yp = cumsum(-sign(y(:,k)));
   y(:,k+1) = 1 + cumsum(yp);
   err = norm(y(:,k+1)-y(:,k)); k = k+1;
end

plot(y(:,end))%, hold on, plot(diff(y(:,end))) % IT WORKS WELL!!!!


%% Now let's try something a little harder: y'' + |y|y = 0
%clear all, close all
dom = [0,6]; tol = 1e-8; err = 2*tol; k = 1;

y = chebfun(@(t)1,dom,'splitting','on');

while err > tol
   yp = 1 + cumsum(-abs(y(:,k)).*y(:,k));
   y(:,k+1) = 1 + cumsum(yp);
   err = norm(y(:,k+1)-y(:,k)); k = k+1;
end

plot(y(:,end))%, hold on, plot(diff(y(:,end)))

%% Repeat with ode45
dom = [0,6]; IC = [1;1];
f = @(t,x) [ x(2); -abs(x(1)).*x(1) ];
[t,x] = ode45(f,dom,IC);
plot(t,x(:,1))

%% Trying to split the domain in two
clear all, close all
dom1 = [0,3]; dom2 = [3,6]; tol = 1e-13; err = 2*tol; k = 1;

y1 = chebfun(@(t)1,dom1,'splitting','on');
y2 = chebfun(@(t)1,dom2,'splitting','on');

while err > tol
   yp1 = 1 + cumsum(-abs(y1(:,k)).*y1(:,k));
   y1(:,k+1) = 1 + cumsum(yp1);
   
   yp2 = yp1(3) + cumsum(-abs(y2(:,k)).*y2(:,k));
   y2(:,k+1) = y1(3,k) + cumsum(yp2);
   
   err = max(norm(y1(:,k+1)-y1(:,k)),norm(y2(:,k+1)-y2(:,k))); k = k+1;   
end

soln = join(y1(:,end),y2(:,end));
plot(soln), ylim([-3,3])
%%
for i = 1:size(y1,2)
   plot(join(y1(:,i),y2(:,i))), ylim([-3,3]), pause(1) 
end

%% Now with AA
clear all, close all
dom = [0,4];
g = @(x) [1 + cumsum(x(:,2)),1 + cumsum(-abs(x(:,1)).*x(:,1))];
init = [chebfun(@(t)1,dom),chebfun(@(t)0,dom)];
soln = AA(g,init);
plot(soln(:,1))


