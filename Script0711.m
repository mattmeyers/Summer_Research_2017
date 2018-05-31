
V = chebop(@(x, u) u + x.*volt(@(x, y) exp(x-y), u), [0, 2]); 
x = chebfun('x',[0,2]);
u = V \ sin(exp(3*x)); 
%%
t = chebfun('t',[0,10]);
V = chebop(@(t,u) 4*u + volt(@(t,s) sin(t-s),u),[0,10]);
u = V \ 5*t;
plot(u)

%%
hold on
plot(t+(1/(2*sqrt(5)))*sin((sqrt(5)/2)*t))

%%
N = chebop(@(t,u)diff(u,2)+u,[0,10]);
N.lbc = [0;5/4];
u = N \ (5/4)*t;
plot(u)

%% Figuring out Volt
clear all, close all
n = 4;
[t,w] = clencurt(n);
A = []
for i = 1:n+1
    c = (t(i)-0)/2;
   for j = 1:n+1
       
   end
end