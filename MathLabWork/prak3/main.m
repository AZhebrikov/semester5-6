c=1;
k=0.04;
m=1;
mu=0.5;
g=9.8;
v=4;
global  m c k mu g sgn v

dt=0.00001;
T0=0;
Tf=50;
X0=[-5,11];
options=odeset('Events','on');

figure(1)
cla
hold on
figure(2)
cla
hold on

sgn=sign(X0(2)-v);
while T0<Tf
    [T,Y,TE,YE,IE]=ode45('bodytr',[T0,Tf],X0,options);

    figure(1)
    plot(T,Y)
    figure(2)
    plot(Y(:,1),Y(:,2))
    
    if (T(end)>=Tf)|(isempty(IE)==1)
        break
    end
    
    T0=T(end);
    X0=Y(end,:);
    q=-c*X0(1)-k*X0(2);
    
    if abs(q)<mu*m*g
        Tr=((1.8*mu*m*g-k*v)/c-X0(1))/v;
        if Tr>0
            t=T0+linspace(0,Tr);
            Y=[X0(1)+v*(t-T0);v*ones(size(t))];
            figure(1)
            plot(t,Y','-.')
            figure(2)
            plot(Y(1,:),Y(2,:),'-.')
            T0=t(end);
            X0=Y(:,end)';
        end
        T0=T0+dt;
        X0=X0+dt*bodytr(t,X0,'')';
        sgn=sign(X0(2)-v);
    end
end
figure(1)
hold off
figure(2)
hold off
