load("no_white.mat");
load("white.mat");
load("trj.txt")
%%
R = 6371e3;
d2r = pi/180;

fsz = 20; font = 'Times'; lw = 2;

subplot(2,3,1);
plot(imu_w(:,1),trj_out_nw(:,8));
grid on;
subplot(2,3,2);
plot(imu_w(:,1),trj_out_nw(:,9));
grid on;
subplot(2,3,3);
plot(imu_w(:,1),trj_out_nw(:,10));
grid on;
subplot(2,3,4);
plot(imu_w(:,1),trj_out_nw(:,8)-trj_out_w(:,8));
grid on;
subplot(2,3,5);
plot(imu_w(:,1),trj_out_nw(:,9)-trj_out_w(:,9));
grid on;
subplot(2,3,6);
plot(imu_w(:,1),trj_out_nw(:,10)-trj_out_w(:,10));
grid on;

clear()