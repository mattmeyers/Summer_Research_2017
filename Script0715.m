%% Let's look into splitting intervals
clear all, close all
dom1 = [0,3]; dom2 = [3,6]; tol = 1e-13; err = 2*tol; k = 1;

y1 = chebfun(@(t)1,dom1,'splitting','on');
y2 = chebfun(@(t)1,dom2,'splitting','on');

while err > tol
   yp1 = 1 + cumsum(-abs(y1(:,k)).*y1(:,k));
   y1(:,k+1) = 1 + cumsum(yp1);
   
   err = norm(y1(:,k+1)-y1(:,k)); k = k+1;   
end

k = 1; err = 2*tol;
while err > tol
   yp2 = yp1(3) + cumsum(-abs(y2(:,k)).*y2(:,k));
   y2(:,k+1) = y1(3,end) + cumsum(yp2);
   
   err = norm(y2(:,k+1)-y2(:,k)); k = k+1;   
end

soln = join(y1(:,end),y2(:,end));
plot(soln), ylim([-3,3])


%% Now let's solve y''+y'+abs(y)y=0
%clear all, close all
dom1 = [0,2]; dom2 = [2,4]; tol = 1e-12; err = 2*tol; k = 1;

y1 = chebfun(@(t)1,dom1,'splitting','on');
y2 = chebfun(@(t)1,dom2,'splitting','on');

while err > tol
   yp1 = 1 + cumsum(-diff(y1(:,k))-abs(y1(:,k)).*y1(:,k));
   y1(:,k+1) = 1 + cumsum(yp);
   err = norm(y1(:,k+1)-y1(:,k)); k = k+1;
end
k=1; err = 2*tol;
while err > tol
   yp2 = yp1(2) + cumsum(-diff(y2(:,k))-abs(y2(:,k)).*y2(:,k));
   y2(:,k+1) = y1(2,end) + cumsum(yp2);
   err = norm(y2(:,k+1)-y2(:,k)); k = k+1;
end
soln = join(y1(:,end),y2(:,end));
plot(soln)