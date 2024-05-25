
%% P

A=[...
    0,0,1,0;...
    0,0,0,1;...
    1/2*9.8-1/8,1/8,0,0;...
    1/8,1/2*9.8-1/8,0,0 ...
    ];

B=[0;0;0;1];
C=[0 1 0 0];
D=0;

sys=ss(A,B,C,D);

U=ctrb(sys);rankU=rank(U);
% the system is fully controllable.
V=obsv(sys);rankV=rank(V);
% the system is fully observable.

P=[-2.3,-2.2,-2.1,-2.0];

K=place(A,B,P);
% u=Kx', where x' - assessment of the system position.
LT = place(A',C',P);
L=LT';
% \dot x' =A x' + B u + L( z - C x' - D u )
% \dot \delta x = (A - L C) \delta x

%% L

%{
 \dot x = A x + B_L u + G w
 \dot z = C x + D_L u + H w + v
    M[w(t)]=0
    M[w(t)w(s)^T] = Q \delta(t-s)
    M[w(t)v(s)^T] = N \delta(t-s)
    M[v(t)]=0
    M[v(t)v(s)^T] = R \delta(t-s)

u=Ko x';
\dot x' = A x' + B_L u + Lo(z − C x' − D_L u)
%}
Q=1;
G=ones(4,1);
R=1;
H=1;
N=0;

sys1=ss(A,[B G],C,[D H]);

[kest,Lo,P]=kalman(sys1,Q,R,N);
[Ko,S,e] = lqr(A,B,Q,R,N);

%% Conclusions

%The control set by the method of setting the eigenvalues of
%P=[-2.3,-2.2,-2.1,-2.0] gave the control matrix K=1.0e+03*[1.4111,0.0373,0.6457,0.0086] and the evaluation matrix L=1.0e+03*[0.6457;0.0086;1.4111;0.0373].
 
%The search for linear quadratic-Gaussian control and the estimation
%algorithm gives: Ko=1.0e+03*[1.4896,0.0389,0.6816,0.0089] and
%Lo=1.0e+03*[0.6388;0.0086;1.3961;0.0370].
   
