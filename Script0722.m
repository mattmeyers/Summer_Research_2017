%% Trying Newton with DDE: y'(t)=-.2y(t-1), phi(t)=1
clear all, close all
dom = [0,5]; hdom = [-1,0];
hist = chebfun(@(t)1,hdom,'splitting','on'); 
yh = chebfun(@(t)1,dom,'splitting','on'); yh = join(hist,yh);

N = chebop(@(t,d)d + (1/5)*volt(@(t,s)1,d(t-1)),[0,5]);
r = 1-cumsum(chebfun(@(t)yh(t-1),dom))-restrict(yh,dom);

y = N\r;
y = yh + y;