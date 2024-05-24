n = input('Введите размер матрицы:');

p=linspace(-20,20,100);
%строю  матрицу корней
res = arrayfun(@(p)(eig_matrix(p,n)),p,'uniform',0);
res = [res{:}]';
roots=res(:,1:end-2)';
conds=res(:,end-1)';
ranks=res(:,end)';
%связываю корни в кривые
res=fun_sort(roots)';


subplot(1,2,1);
plot(res); 
grid on;
xlabel('Re $\lambda$','Interpreter','latex');
ylabel('Im $\lambda$','Interpreter','latex');

subplot(1,2,2);
plot(p',conds); 
grid on;
xlabel('parametr, $p$','Interpreter','latex');
ylabel('$cond(p)$','Interpreter','latex');

clear()

%%
function y = eig_matrix(p,n)
A = construction_matrix(p,n);
y=eig(A);
y = [y;cond(A);rank(A)];
end
% преобразовал уравнение к линейному и работаю с полученной матрицей.
function A = construction_matrix(p,n)
    A=eye(2*n);
    A(1:n,1:n)=zeros(n);
    A(1:n,(n+1):2*n)=-eye(n);
    A((n+1):2*n,1:n)=(2*eye(n)+diag(ones(1,n-1),1)+diag(ones(1,n-1),-1))+(diag((-1).^(1:1:n-1),1)+diag((-1).^((1:1:n-1)+1),-1)).*p;
    A((n+1):2*n,(n+1):2*n)=0.1*(2*eye(n)+diag(ones(1,n-1),1)+diag(ones(1,n-1),-1));
end
% функция поиска собственных корней возвращает корни не в установившемся
% порядке, поэтому произвожу следующую сортировку:
% считая вектор корней при предыдущем параметре верно упорядоченным,
% рассматриваю матрицу всех возможных перестановок последующего вектора
% корней, и выбираю такой вектор, которой менее всех отстоит от вектора при прошлом параметре. 
function A = fun_sort(B)
A(:,1)=B(:,1);
    for k=2:size(B,2)
        tmp_matrix=(perms((B(:,k))'))'-A(:,k-1);
        min_vector=tmp_matrix(:,1);
        for s=2:size(tmp_matrix,2)
            if norm(tmp_matrix(:,s),1)<norm(min_vector,1)
                min_vector=tmp_matrix(:,s);
            end
        end
        A(:,k)=A(:,k-1)+min_vector;
    end
end
