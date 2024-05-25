sinit.x=linspace(0,pi,200);
sinit.y=[-1+2/pi*sinit.x;zeros(size(sinit.x))];

f=@(t,u)[u(2);-u(1)];
g=@(ya,yb)[ya(1)+1;yb(1)-1];
sol=bvp4c(f,g,sinit);
plot(sol.x,sol.y(1,:));

[X,Y] = ode45(f,[0,pi],sol.y(:,1),options);
Y=sol.y';