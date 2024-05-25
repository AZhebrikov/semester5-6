

%{
M=linspace(0,1,1000);
M=M';
[K,E]=ellipke(M);
y=2.*E./K-1;
E_d=1/2./M.*(E-K);
K_d=1/2./M./(1-M).*E-1/2./M.*K;
F=4.*(y.^3)./3.*( (K.^2)./y - (K.^3).*K_d./(E_d.*K-K_d.*E));
plot(M,F,LineWidth=3)
ylim([-5,10])
grid on
%}


L=2;
I=1;
E=1;
G=6.67430e-11;
alpha=1.01;
m=sqrt(alpha*E*I*pi^2/G);

y=linspace(0.001,L/2,1000);
W=y;
P=y;
F_w=y;
F=y;
P_y=y;
M=y;
k=y;
V_yy=y;

for l=1:length(y)

y0 = y(l);
M0 = fun_y(y0,L);
M(l)=M0;

[K0,E0] = ellipke(M0);
P0=4*I*E/(L^2)*K0^2;
P(l)=P0;

F0=G*m^2/4/y0^2;
F(l)=F0;

w0=G*m/(4*y0^3)-P0/(m*y0);
W(l)=w0;

F_w0=m*y0*w0;
F_w(l)=F_w0;

P_y0=4*P0/L*(E0-(1-M0)*K0)/((1-2*M0)*E0-(1-M0)*K0);
P_y(l)=P_y0;

V_yy0=-3/4*G*m^2/y0^3+P0/y0-P_y0;
V_yy(l)=V_yy0;

end    

figure()
fig=nexttile;
hold on
%plot(y',W',LineWidth=2,Color='#0072BD')
plot(y',V_yy',LineWidth=2,Color='#A2142F')
%plot(y',V_yy',LineWidth=2,Color='#77AC30')
hold off
lgd=legend({'$V_{yy}$'},Interpreter="latex")
ylim([-100,10])
grid on

function fun_out = fun_y_M(M,L)
   [K,E] = ellipke(M);
   fun_out = L./2.*(2.*E./K-1);
end

function fun_out = fun_y(y0,L)
    M1=0.5;
    M2=0.5;
        
    while fun_y_M(M1,L)-y0<0
        M1=M1/2;
    end
    while fun_y_M(M2,L)-y0>0
        M2=(M2+1)/2;
    end

    while M2-M1>10e-9
        if fun_y_M((M1+M2)/2,L)-y0 >0
            M1=(M1+M2)/2;
        elseif fun_y_M((M1+M2)/2,L)-y0 <0
            M2=(M1+M2)/2;
        elseif fun_y_M((M1+M2)/2,L)-y0 == 0
            break
        end
    end
    fun_out=(M1+M2)/2;
end




