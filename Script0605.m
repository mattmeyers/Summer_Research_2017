%% Implementing Picard iteration to find convergence of y' = L*y
clear all, close all
L = 8;  y0 = 1;  tspan = [0,5];

%% Finding exact solution
M = chebop(tspan);
M.op = @(t,y) diff(y)-L*y;
M.lbc = y0;
yExact = M\0;
plot(yExact)

%% Performing Picard iteration
y = chebfun('y',tspan);
numiter = 206;

for i = 1:numiter
    yApp = y0 + cumsum(L*y);
    y = yApp;
end

plot(yApp)

%% Comparing solutions
%clf
%plot(yExact)
%str = sprintf('Exact solution: L=%d, y0=%d',L,y0);
%title(str)
%print -depsc L8Exact.eps
clf
str = sprintf('Approximate solution: L=%d, y0=%d, numiter=%d',L,y0,numiter);
plot(yApp), title(str)
print -depsc L8App206.eps


