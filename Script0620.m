%% Repeating DDEs with Anderson Acceleration [y'(t)=-y(t-1)]
clear all, close all
dom = [0,4];
hist = chebfun(@(t)(-1).^floor(-4*t),[-1,0],'splitting','on');
x = chebfun(@(t) 1,[0,4]);
y = join(hist,x);

g = @(x) hist(0) + cumsum(-y(x-1),dom);

for i = 1:10
    x = g(x);
end

plot(x)