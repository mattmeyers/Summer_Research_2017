%% Harmonic Oscillator
clear all, close all
dom = [0,5];
g = @(x) [1 + cumsum(chebfun(@(t) x(t,2),dom,'splitting','on')),cumsum(chebfun(@(t) -x(t,1),dom,'splitting','on'))];
x = [chebfun(@(t) 1, dom,'splitting','on'), chebfun(@(t) 0, dom,'splitting','on')];

[soln,x,numIter,errVec] = AA(g,x);
plot(soln)


%% A more interesting problem: y''+sign(y)=0
dom = [0,5];
N = chebop(dom);
N.op = @(t,y) diff(y,2)+sign(y);
%N.lbc = 1;
y = chebfun.ode113(N,domain(dom),1);
%y = N\0;  % Doesn't work...
plot(y)

%% Solving with AA now
clear all, close all
dom = [0,5];
g = @(x) [1 + cumsum(x(:,2)),cumsum(-sign(x(:,1)))];
init = [chebfun(@(t)1,dom),chebfun(@(t)0,dom)];
[soln,x,numIter,errVec] = AA(g,init);
plot(soln)

%% What about with ode45
clear all
dom = [0,20]; IC = [1,0];
f = @(t,x) [x(2); -sign(x(1))];

[t,x] = ode45(f,dom,IC);
plot(t,x)


%% Moving onto a first order equation, first with ode45
clear all, close all
dom = [0,5];
f = @(t,x) -abs(x)*x + 3*cos(t);

[t,x] = ode15s(f,dom,0); 
plot(t,x)

%% Time to try Chebfun
clear all, close all
dom = [0,5];
N = chebop(@(t,y)diff(y)+sign(y).*y, dom);
N.lbc = 0; rhs = chebfun('3*cos(t)',dom);
y = N\rhs;
plot(y)

%% And now AA
clear all, close all
dom = [0,5];
g = @(x) cumsum(chebfun(@(t) -abs(x(t)).*x(t)+3*cos(t),dom,'splitting','on'));
init = chebfun(@(t)sin(t),dom);

[soln,~,numIter] = AA(g,init);
plot(soln,'-o')

