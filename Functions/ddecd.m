function [y] = ddecd(ddefun,dom,hist,init)
a = dom(1); b = dom(2);
ya = hist(a);

y = join(hist,init);

err = 1; tol = 1e-8;

while err > tol
   yNew = ya + cumsum(chebfun(@(t)ddefun(t,y),dom,'splitting','on'));
   yNew = join(hist,yNew);
   err = norm(y-yNew);
   y = yNew;
end

end