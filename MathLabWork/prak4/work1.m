% x'=Ax + Bu
% z = Cx + Du
A = [1, eps; 0, 1];
B = [1; eps];
rank(ctrb(A, B));
ctrb(A, B);
obsv(A, B)