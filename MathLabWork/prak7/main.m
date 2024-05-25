M=0.04;
m=0.312;
Jm=0.0000803;
c1=0.0069;
c2=0.000099;
u0=19;
g=9.8;
l=0.3;
J=M*l^2/2;
sigma_mom=0.3;

a1=(M/2+m)/(J+m*l^2)*g*l;
a2=(J+m*l^2+Jm)/(J+m*l^2)/Jm;
c11=c1/(J+m*l^2);
c21=c2/(J+m*l^2);
c12=a2*c1;
c22=a2*c2;
A=[0 1 0 0; a1 0 0 c21; 0 0 0 1; -a1 0 0 -c22];
B=[0 0; 1/(J+m*l^2) -c11 ; 0 0;  0 -c12];
C=[0 0 0 1;1 0 0 0;0 0 1 0];
D=[0 0; 0 0; 0 0];

sys=ss(A,B,C,D);
sysop=sys(:,1)*tf(1,[0.1,0]);
U=ctrb(sysop);
rank(U)
v=obsv(sysop(2,:));
rank(v)
K=place(sysop(1,1),[-1,-1.2,-5,-5.7,-10]);
L=place(syaop.a',sysop.c(2,:)',[-1,-1.2,-5,-5.7,-10]);
sysr=ss(sysop.a-sysop.b(:,1)*K-L'*sysop.c(2,:),L',-K,0);
syscl=feedback(sysop,sysr,1,2,1);
eig(syscl,a)

figure(1)
[Y,T,X]=initial(syscl,[1 0 0 0 0 0 0 0 0 0]);
plot(T,X(:,1)-X(:,6))