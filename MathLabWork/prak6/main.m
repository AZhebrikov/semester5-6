global u


dt=0.00001;
eps=1e-9;
Eps=1;
options = odeset('Events','on','RelTol',1e-9,'AbsTol',1e-12);

figure(1)
cla
hold on
figure(2)
cla
hold on

T0=0;
X0=[3,1];
while 1>0

    control_selection(T0,X0);
    
    if norm(X0,T0)<eps
        [T,Y,TE,YE,IE]=ode45('fun_right_eps',[T0,T0+100],X0,options);
    else
        [T,Y,TE,YE,IE]=ode45('fun_right',[T0,T0+100],X0,options);
    end

    figure(1)
    plot(T,Y)
    grid on
    figure(2)
    plot(Y(:,1),Y(:,2))
    grid on

    T0=TE;
    X0=YE;

    if max(abs(X0(1)),abs(X0(2)))<Eps
        break
    end
end

T0=0;
X0=[3,1];
while 1>0

    control_selection(T0,X0);
    
    [T,Y,TE,YE,IE]=ode45('fun_right_v',[T0,T0+100],X0,options);
    
    figure(1)
    plot(T,Y)
    grid on
    figure(2)
    plot(Y(:,1),Y(:,2))
    grid on

    T0=TE;
    X0=YE;

    if max(abs(X0(1)),abs(X0(2)))<Eps
        break
    end
end

figure(1)
hold off
figure(2)
hold off

function f_out = f_switching_line_x(t,x)
    if x(2)>=0
        f_out=x(1)+x(2)^2/2;
    else
        f_out=x(1)-x(2)^2/4;
    end
end

function f_out = f_switching_line_y(t,x)
    if x(1)>=0
        f_out=x(2)+sqrt(4*x(1));
    else
        f_out=x(2)-sqrt(-1*2*x(1));
    end
end


function u_out = control_selection(t,x)
global u
    if f_switching_line_x(t,x) > 0
        u=-1;
    elseif f_switching_line_x(t,x) < 0
        u=2;
    else
        if x(2)>0
            u = -1;
        else 
            u = 2;
        end
    end
u_out=0;
end

function norm_out = norm(x,t)
norm_out = max(abs(f_switching_line_x(t,x)),abs(f_switching_line_y(t,x)) );
end