%% Figuring out how to determine when to split
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
%% Looking at the norms
norms = zeros(15,2);
for i = 1:15
   norms(i,1) = norm(y1(:,i));
   norms(i,2) = norm(y2(:,i)); %norm gets larger then NaN
end

%% Trying to implement auto splitting
dom = [0,4]; ic = [1;1];
g = @(y) -abs(y).*y;
y = chebfun(@(t)1,dom);

tol = 1e-8; err = 2*tol; k = 1;
while err > tol
   z = ic(1) + cumsum(g(y(:,k)));
   for i = 2:length(ic)
       z(:,i) = ic(i) + cumsum(z(:,i-1));
   end
   y(:,k+1) = z(:,end);
   err = norm(y(:,k+1)-y(:,k)); k = k+1;
end
plot(y(:,end))

%% Let's try the new function
clear all, close all
dom = [0,12]; ic = [1;1];
g = @(y) -abs(y).*y;
y = chebfun(@(t)1,dom);

soln = intsplit(g,y,ic);
plot(soln)

%% Another example
clear all, close all
dom = [0,5]; ic = [0;1];
g = @(y)-sign(y).*y;
y = chebfun(@(t)cos(t),dom,'splitting','on');

soln = intsplit(g,y,ic);
plot(soln)

%% Yet another
clear all, close all
dom = [0,20]; ic = 0;
t = chebfun('t',dom);
g = @(y)-y.^3+3*cos(t);
y = chebfun(@(t)cos(t),dom);

soln = intsplit(g,y,ic);
plot(soln)

%% Kermack McKendrick, intervals
clear all, close all
dom = [0,1]; IC = [5,.1]; R = 1;

x1 = chebfun(5,dom);
x2 = chebfun(.1,dom);
x3 = chebfun(1,dom);
init = [x1,x2];

g = @(x) [IC(1) + cumsum(-R*x(:,1).*x(:,2)), ...
    IC(2) + cumsum(R*x(:,1).*x(:,2)-x(:,2))];

soln = AA(g,init);




