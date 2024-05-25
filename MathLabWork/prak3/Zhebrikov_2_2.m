solinit.x=linspace(0,10,1000);
solinit.y=[zeros(size(solinit.x));ones(size(solinit.x))];

f = @(t,u)[u(2);-0.1*u(2)-(1+0.1*sin(2*t))*u(1)];
bc=@(x_init,x_end)[x_init(1);x_end(1)-x_end(2)-4];
sol=bvp4c(f,bc,solinit);

t=sol.x';
x=sol.y';
bc=[x(:,1),x(:,1)-x(:,2)-4];

options=odeset('AbsTol',1e-15,'RelTol',1e-9);
[T,X] = ode45(f,[0,10],x(1,:),options);
BC=[X(:,1),X(:,1)-X(:,2)-4];

t_inter=linspace(0,10,1000);
x_inter = interp1(t',x(:,1)',t_inter,'makima');
X_inter = interp1(T',X(:,1)',t_inter,'makima');
t_inter=t_inter';
x_inter=x_inter';
X_inter=X_inter';

subplot(1,3,1)
grid on
hold on
plot(t,x)
plot(T,X,'-.')
legend({'$x(t) (bvp4c)$','$\dot x(t) (bvp4c)$','$x(t) (ode45)$','$\dot x(t) (ode45)$'},'Interpreter','latex')
xlabel('t')

subplot(1,3,2)
grid on
plot(t_inter,x_inter-X_inter)
xlabel('t')
ylabel('$\Delta x = x(bvp4c) - x(ode45)$','Interpreter','latex')

subplot(1,3,3)
grid on
hold on
plot(t,bc)
plot(T,BC,'-.')
legend({'$bc1(t) (bvp4c)$','$bc2(t) (bvp4c)$','$BC1(t) (ode45)$','$BC2(t) (ode45)$'},'Interpreter','latex')
xlabel('t')
