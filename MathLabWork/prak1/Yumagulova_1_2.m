clear;
M = 1;
c1 = linspace(1,3,30);
n = 5;
c = 1;
k = 0.001;

roots=zeros(2*n,length(c1));
for i = 1:length(c1)

    B=diag(c*ones(n-1,1),-1)+diag(c*ones(n-1,1),1)+diag(-2*c*ones(n,1));
    B(1,1) = -c1(i)-c;
    B=B./M;

    C=diag(-1*k*ones(n,1));
    C=C./M;

    A = [zeros(n),eye(n);C,B];
    roots(:,i) = eig(A);
end

roots=roots*(-sqrt(-1));
roots=sort(roots,'ComparisonMethod','real');
roots=roots*(sqrt(-1));

roots=roots';
figure(1);
clf;
plot(c1',roots);

xlabel('Re');
ylabel('Im');
title('Roots trajectories with changing hardness c1');
%legend(cellstr(num2str(c1)));
grid on;
