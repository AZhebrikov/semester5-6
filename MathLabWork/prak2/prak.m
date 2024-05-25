
load("imu.txt")
load("trj.txt")
%% постоянные

a=6378137;
e_2=6.6943799901413*10^(-3);
u = 7.292115e-5;
dt=0.01;
d2r = pi/180;

%% начальная выставка

% время неподвижного пребывания где-то в районе 90 секунд, по которым и
% проведу измерение.(9000 итераций записи).

matrix_initial_exhibition=imu(1:9000,2:end);
v_initial_exhibition=mean(matrix_initial_exhibition);
v_in_angular_velocity = d2r*v_initial_exhibition(1:3)';
v_in_specific_force = v_initial_exhibition(4:6)';

L_initial=[...
    cross(v_in_angular_velocity,v_in_specific_force)/norm(cross(v_in_angular_velocity,v_in_specific_force)),...
    cross(v_in_specific_force,cross(v_in_angular_velocity,v_in_specific_force))/norm(cross(v_in_specific_force,cross(v_in_angular_velocity,v_in_specific_force))),...
    v_in_specific_force/norm(v_in_specific_force)];

%% начальные условия

fi_now = trj(9000,2)*d2r;
lambda_now = trj(9000,3)*d2r;
h_now = trj(9000,4);

V_now = zeros(3,1);

psi_in = atan2(L_initial(1,1),L_initial(1,2));
tetha_in = atan2( L_initial(1,3) , sqrt( (L_initial(1,1))^2+(L_initial(1,2))^2 ) );
gamma_in = atan2( -L_initial(3,3) , L_initial(2,3) );

    w_vec_now = v_in_angular_velocity;
    f_vec_now = v_in_specific_force;

        A_x_matrix_now = eye(3); 
        A_z_matrix_now = L_initial;


trj_out = trj;

trj_out(9000,2) = fi_now/d2r; 
trj_out(9000,3) = lambda_now/d2r; 
trj_out(9000,4) = h_now;

trj_out(9000,5:7) = V_now'; 

trj_out(9000,8) = psi_in/d2r; 
trj_out(9000,9) = tetha_in/d2r; 
trj_out(9000,10) = gamma_in/d2r;    

%% интегрирование
for t=9000:1:(size(imu,1)-9000)
    
    R_N = a*(1-e_2)/(sqrt( 1-e_2*(sin(fi_now))^2 )^3);
    R_E = a/sqrt( 1-e_2*(sin(fi_now))^2 );
    
    OMEGA = [-V_now(2)/(R_N+h_now),V_now(1)/(R_E+h_now),V_now(1)*tan(fi_now)/(R_E+h_now)];
    OMEGA_ = [0,OMEGA(3),-OMEGA(2);-OMEGA(3),0,OMEGA(1);OMEGA(2),-OMEGA(1),0];
    
    U = [0,u*cos(fi_now),u*sin(fi_now)];
    U_ = [0,U(3),-U(2);-U(3),0,U(1);U(2),-U(1),0];
    
    W_x = OMEGA_ + U_;
    A_x_matrix_next = (  eye(3)+sin(norm(W_x)*dt)/norm(W_x)*W_x+(1-cos(norm(W_x)*dt))/((norm(W_x))^2)*W_x*W_x  )*A_x_matrix_now;   
    
    W_z = [0,w_vec_now(3),-w_vec_now(2);-w_vec_now(3),0,w_vec_now(1);w_vec_now(2),-w_vec_now(1),0];
    A_z_matrix_next = (  eye(3)+sin(norm(W_z)*dt)/norm(W_z)*W_z+(1-cos(norm(W_z)*dt))/((norm(W_z))^2)*W_z*W_z  )*A_z_matrix_now;   
    
    L_next = A_z_matrix_next*(A_x_matrix_next)';
    psi_next = atan2(L_next(1,1),L_next(1,2));
    tetha_next = atan2( L_next(1,3) , sqrt( (L_next(1,1))^2+(L_next(1,2))^2 ) );
    gamma_next = atan2( -L_next(3,3) , L_next(2,3) );
    
    fi_next = fi_now + V_now(2)/(R_N+h_now)*dt;
    lambda_next = lambda_now + V_now(1)/((R_E+h_now)*cos(fi_now))*dt;
    h_next = h_now + V_now(3)*dt;
    
    L_now = A_z_matrix_now*(A_x_matrix_now)';
    V_next = V_now + dt.*( (OMEGA_ + 2.*U_)*V_now + g_(fi_now,h_now) + L_now'*f_vec_now );
    
    w_vec_next = d2r*imu(t+1,2:4)';
    f_vec_next = imu(t+1,5:7)';

    %-----

    trj_out(t+1,2) = fi_next/d2r; 
    trj_out(t+1,3) = lambda_next/d2r; 
    trj_out(t+1,4) = h_next;

    trj_out(t+1,5:7) = V_next'; 
    
    trj_out(t+1,8) = psi_next/d2r; 
    trj_out(t+1,9) = tetha_next/d2r; 
    trj_out(t+1,10) = gamma_next/d2r;

    %---

    A_x_matrix_now = A_x_matrix_next; 
    A_z_matrix_now = A_z_matrix_next;
    V_now = V_next;
    fi_now = fi_next;
    lambda_now = lambda_next;
    h_now = h_next;

    w_vec_now = w_vec_next;
    f_vec_now = f_vec_next;

end
%% вывод данных
imu_nw=imu;
trj_out_nw=trj_out;
save("no_white.mat","imu_nw","trj_out_nw",'-mat');

%% вывод данных
fsz = 20; font = 'Times'; lw = 2;
trj(:,2:3) = trj(:,2:3)*d2r;
trj_out(:,2:3) = trj_out(:,2:3)*d2r;

R = 6371e3;

%plot(imu_save(:,1),imu_save(:,2:4)-imu(:,2:4));
%grid on;

subplot(3,2,1);

hold on;
plot(trj(1:end,1),[trj_out(1:end,2)*R,trj(1:end,2)*R],'LineWidth',lw);
plot(trj(1:end,1),[trj_out(1:end,3)*R.*cos(trj(1,2)),trj(1:end,3)*R.*cos(trj(1,2))],'LineWidth',lw); 
plot(trj(1:end,1),[trj_out(1:end,4),trj(1:end,4)],'LineWidth',lw);
grid on;
xlabel('Время, с'); ylabel('Координаты, м');
legend('x_1','x_2','x_3');
set(gca,'FontSize',fsz,'FontName',font,'LineWidth',lw);

subplot(3,2,2);

hold on;
plot(trj(1:end,1),(trj_out(1:end,2)-trj(1:end,2))*R,'LineWidth',lw);
plot(trj(1:end,1),(trj_out(1:end,3)-trj(1:end,3))*R.*cos(trj(1,2)),'LineWidth',lw); 
plot(trj(1:end,1),trj_out(1:end,4)-trj(1:end,4),'LineWidth',lw);
grid on;
xlabel('Время, с'); ylabel('Ошибка координаты, м');
legend('\Delta x_1','\Delta x_2','\Delta x_3');
set(gca,'FontSize',fsz,'FontName',font,'LineWidth',lw);

subplot(3,2,3);

plot(trj(1:end,1), [trj_out(1:end,5:7),trj(1:end,5:7)],'LineWidth',lw);
legend('v_1','v_2','v_3');
xlabel('Время, с'); ylabel('Cкорости, м/с');
set(gca,'FontSize',fsz,'FontName',font,'LineWidth',lw);

subplot(3,2,4);

plot(trj(1:end,1), trj_out(1:end,5:7)-trj(1:end,5:7),'LineWidth',lw);
legend('\Deltav_1','\Deltav_2','\Deltav_3');
xlabel('Время, с'); ylabel('Ошибка скорости, м/с');
set(gca,'FontSize',fsz,'FontName',font,'LineWidth',lw);

subplot(3,2,5);

plot(trj_out(:,1),[trj_out(:,8:10),trj(:,8:10)]);
legend('\phi','\theta','\gamma');
xlabel('Время, с'); ylabel('Углы, м/с');
set(gca,'FontSize',fsz,'FontName',font,'LineWidth',lw);

subplot(3,2,6);

plot(trj(:,1),trj_out(:,8:10)-trj(:,8:10));
legend({'\Delta \phi','\Delta \theta','\Delta \gamma'});
xlabel('Время, с'); ylabel('Ошибка угла, м/с');
set(gca,'FontSize',fsz,'FontName',font,'LineWidth',lw);

clear()
%% функции
function g = g_(fi,h)
a = 6378137;
g_e = 9.78030;
betta_1 = 5.302*10^(-3);
betta_2 = 7*10^(-6);
delta_d = 14*10^(-5);
g=[0,0,-g_e*( 1+betta_1*(sin(fi))^2-betta_2*(sin(2*fi))^2 - 2*h/a ) + delta_d ]';
end