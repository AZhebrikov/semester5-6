clearvars;
%% configuration
imu_file = 'imu.txt';
trj_file = 'trj.txt';
tc = 1; omc = 2:4; fc = 5:7;
llhc = 2:4; vc = 5:7; eac = 8:10; 
fs = 100; % sampling frequency, Hz
t0 = 90; % alignment time, seconds
u =  7.292115e-5; % angular rate of Earth, rad/s
a = 6378137; e2 = 6.6943799901413e-3; % parameters of Earth ellipsoid
d2r = pi/180; % degree to rad
%% loading data
imu = load(imu_file);
t = imu(:,tc); om = imu(:,omc)*d2r; f = imu(:,fc);
N = size(imu,1); % number of samples
trj = load(trj_file);
llh_trj = trj(:,llhc); % trajectory latitude, longitude, height
v_trj = trj(:,vc); % trajectory velocity
ypr_trj = trj(:,eac); % trajectory yaw, pitch, roll
%% initial alignment
L0 = initial_alignment(t0,om,f,fs);
ypr0 = get_ypr(L0);
%% main navigation cycle
% initializing structures for navigation solution
llh = zeros(N,3); v = zeros(N,3); ypr = zeros(N,3);
llh(1:t0*fs+1,:) = llh_trj(1:t0*fs+1,:);
llh(:,1:2) = llh(:,1:2)*d2r;
ypr(1:t0*fs+1,:) = repmat(ypr0,t0*fs+1,1);
Az0 = L0; Ax0 = eye(3); L = L0;
% determining position, velocity, orientation
for step = 2:N
    % skip alignment
    if((step-1)/fs<=t0)
        continue
    end
    fi = llh(step-1,1); lam = llh(step-1,2); h = llh(step-1,3);
    v1 = v(step-1,1); v2 = v(step-1,2); v3 = v(step-1,3);
    % integrate dynamic equations
    %   right side of dynamic equations
    Re = a/sqrt(1-e2*sin(fi)^2);
    Rn = a*(1-e2)/(sqrt(1-e2*sin(fi)^2))^3;
    Omx = [-v2/(Rn+h),v1/(Re+h),v1*tan(fi)/(Re+h)];
    ux = [0,u*cos(fi),u*sin(fi)];
    f_llh = [v2/(Rn+h), v1/(Re+h)/cos(fi), v3];
    f_v = v(step-1,:)*(skew_symm(Omx)+2*skew_symm(ux))';
    f_v = f_v + [0,0,-helmert(fi,h)];
    f_v = f_v + f(step-1,:)*L;
    %   Euler integration step   
    sol = euler_integrate([llh(step-1,:),v(step-1,:)],[f_llh,f_v],1/fs);
    llh(step,:) = sol(1:3);
    v(step,:) = sol(4:6);
    % integrate Poisson equation
    %   for instrumental frame
    Az = poisson(Az0,om(step-1,:),1/fs);
    %   for geographical frame
    Ax = poisson(Ax0,Omx+ux,1/fs);
    %   compute L matrix
    L = Az*Ax';
    ypr(step,:) = get_ypr(L);
    Az0 = Az; Ax0 = Ax;
end
%% plotting 
fsz = 20; font = 'Times'; lw = 2;
llh_trj(:,1:2) = llh_trj(:,1:2)*d2r;
R = 6371e3;
subplot(1,2,1);
% t = t(1:end-1);
% llh_trj = llh_trj(1:end-1,:);
% llh = llh(2:end,:);
plot(t,(llh(:,1)-llh_trj(:,1))*R,'LineWidth',lw); grid on; hold on;
plot(t,(llh(:,2)-llh_trj(:,2))*R.*cos(llh(1,1)),'LineWidth',lw); 
plot(t,llh(:,3)-llh_trj(:,3),'LineWidth',lw);
xlabel('Time, s'); ylabel('Position error, m');
legend('\Delta x_E','\Delta x_N','\Delta x_U');
set(gca,'FontSize',fsz,'FontName',font,'LineWidth',lw);
subplot(1,2,2);
t = t(1:end-1);
v = v(2:end,:); v_trj = v_trj(1:end-1,:);
plot(t,v-v_trj,'LineWidth',lw); grid on;
ylim([-2 2]);
legend('\DeltaV_E','\DeltaV_N','\DeltaV_U');
xlabel('Time, s'); ylabel('Velocity error, m/s');
set(gca,'FontSize',fsz,'FontName',font,'LineWidth',lw);
%% functions

% compute initial orientation matrix
function L0 = initial_alignment(t0,om,f,fs)
    om0 = mean(om(1:t0*fs,:));
    f0 = mean(f(1:t0*fs,:));
    l1 = cross(om0,f0)/norm(cross(om0,f0));
    l2 = cross(f0,cross(om0,f0))/norm(cross(f0,cross(om0,f0)));
    l3 = f0/norm(f0);
    L0 = [l1',l2',l3'];
end

% integrate Poisson equation from t to t+dt
function L = poisson(L0,om,dt)
    L = eye(3);
    L = L + sin(norm(om)*dt)/norm(om)*skew_symm(om);
    L = L + (1-cos(norm(om)*dt))/norm(om)^2*skew_symm(om)^2;
    L = L*L0;
end

% get yaw, pitch and roll from orientation matrix
function ypr = get_ypr(L)
    y = atan2(L(1,1),L(1,2));
    p = atan2(L(1,3),sqrt(L(2,3)^2+L(3,3)^2));
    r = -atan2(L(3,3),L(2,3));
    ypr = [y,p,r];
end

% create skew-symmetric matrix from vector x
function x_hat = skew_symm(x)
    x_hat = [0 x(3) -x(2);-x(3) 0 x(1);x(2) -x(1) 0];
end

% compute gravity with Helmert formula
function g = helmert(fi,h)
    ge = 9.78030; a = 6378137;
    beta1 = 5.302e-3; beta2 = 7e-6;
    dg = 14e-5;
    g = ge*(1+beta1*sin(fi)^2-beta2*sin(2*fi)^2-2*h/a)-dg;
end

% integrate step of differential equation by Euler method
function x = euler_integrate(x0,f,dt)
    x = x0+f*dt;
end