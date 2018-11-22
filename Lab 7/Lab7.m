%% Simon Honigmann
% Sensor Orientation
% Lab 7: Kalman Filtering with simulated GPS data: simple (a=0) model
% 11/16/18

%% Cleanup
clc;
clear;
close all;

%% Lab Formatting:
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultLineLineWidth',1.5);

%% Constants
r = 25; % Circle radius: 500 m, Angular speed ? = ?/100;
omega = pi/100; % angular speed = pi/100
g = 9.81; %gravity, m/s^2
%only for 100Hz and Trapezoidal
sampling_rate = .5; %is this what is meant bz time interval?
method = 'Trapezoidal';
f = sampling_rate;
sample_rate_text = num2str(sampling_rate,3);
num_rotations = 1; %number of time the vehicle goes around the circle. increase to better show divergence

%% Assumptions/Constraints
dr = 0; %constant radius of motion in m frame
ddr = 0;
dpsi = omega; %constant angular velocity in m frame
ddpsi = 0;

%% Initial Conditions
psi_0 = 0; % Initial position: on North axis

alpha_0 = psi_0 + pi/2;
x_0 = [r,0];
dx_0 = [0,omega*r]; % Initial velocity: north-axis: 0, east-axis: ? ? radius
required_time = 2*num_rotations*pi/omega; %time for 1 full rotation [s]
t = (0:1/sampling_rate:required_time)'; %time vector, [s]
samples = (1:length(t))';

%% Generating Actual Trajectories
psi = (psi_0:num_rotations*2*pi/(length(t)-1):num_rotations*2*pi)'; %polar angle, [rad]
alpha = psi + pi/2; %azimuth angle, Initial azimuth: 90? (towards x-accelerometer)
x_m = [r*cos(psi),r*sin(psi)]; % map frame, [north, east]
dx_m = [dr.*cos(psi)-r.*sin(psi).*dpsi , r.*cos(psi).*dpsi + dr.*sin(psi)];
ddx_m = [ddr.*cos(psi) - dr.*sin(psi).*dpsi - (dr.*sin(psi).*dpsi + r.*cos(psi).*dpsi.^2 + r.*sin(psi).*ddpsi), ... %ddx_m1
    dr.*cos(psi).*dpsi - r.*sin(psi).*dpsi.^2 + r.*cos(psi).*ddpsi + ddr.*sin(psi)+dr.*cos(psi).*dpsi ]; %ddx_m2
dalpha = dpsi*ones(length(t),1); %deriving alpha wrt time = deriving psi wrt time
v_b = [ones(length(t),1)*omega*r,zeros(length(t),1)]; % got lazy... this works here because omega and r are constant. Would need to change this line if that was not the case
a_b = [zeros(length(t),1),ones(length(t),1)*omega^2*r]; %got lazy... this works here because omega and r are constant. Would need to change this line if that was not the case

%% Simulating Noisy GPS Signal
sigma_n = 0.5; % random noise of gps measurement in x direction
sigma_e = 0.5; % random noise of gps measurement in x direction

seed = 1234567;
n = length(t);
[n_gpsN,~] = whiteNoiseRandomWalk(n,seed); %generate white noise process
[n_gpsE,~] = whiteNoiseRandomWalk(n,seed+1); %generate white noise process

n_gps = [sigma_n*n_gpsN,sigma_e*n_gpsE];
gps = x_m+n_gps;

%% Kalman Filter Setup
dt = 2; % time step used in prediction

%initialize x, z, F, H, and G
x = zeros(4,n); % pn, pe, vn, ve
x(:,1) = [x_m(1,:)';dx_m(1,:)'];% assume initial state is known

xtild = 0*x;%assume no initial state knowledge
%xhat = 0*x;%assume no initial state knowledge
xhat = x;%assume perfect initial state knowledge
z = gps'; %measurements, with white noise uncertainty

syms deltat;
Adt = [0 0 deltat 0; 0 0 0 deltat; 0 0 0 0; 0 0 0 0];
F=subs(expm(Adt),deltat,dt);%dynamic matrix calculated using matrix exponential

H = [1 0 0 0; 0 1 0 0]; %measurement design matrix

%noise shaping matrix
%G = [dt*eye(2);eye(2)]; %not actually used here since Qk is pre-derived

%Create R
sigma_x = 0.5;
R = sigma_x*eye(2); %covariance of measurement noise

% Solve for Qk using derivation from lecture
dq_v = 0.05; %sigma_vdot
dq_vn = dq_v;
dq_ve = dq_v;

Qk(1,:,:) = [1/3*dt^3*dq_vn, 0 ,1/2*dt^2*dq_vn, 0;...
      0, 1/3*dt^3*dq_ve, 0 ,1/2*dt^2*dq_ve; ...
      1/2*dt^2*dq_vn, 0,dt*dq_vn, 0; ...
      0, 1/2*dt^2*dq_ve, 0, dt*dq_ve]; %covariance matrix of the system noise

%Define P - uncertainty of inital measurement 
sigma_x_state = 10; 
sigma_v_state = 0.1; 
P(1,:,:) = diag([sigma_x_state,sigma_x_state,sigma_v_state,sigma_v_state]).^2; %covariance of state estimate
%Phat = P; % a posteriori covariance of state estimate (after measurement)
Ptild = P; % a priori covariance of state estimate (before measurement)

I = eye(size(squeeze(P(1,:,:))));
  
%% Kalman Filter Loop
%try different prediction deltaT and update deltaT
k=1;
while k<n
%for k=1:n-1
   
    %prediction:
    xtild(:,k+1)=F*xhat(:,k);
    Ptild(k+1,:,:) = F*squeeze(P(k,:,:))*F'+squeeze(Qk(1,:,:));
    
    k=k+1;
    
    %calculate kalman gain
    K(k,:,:) = squeeze(Ptild(k,:,:))*H'/(H*squeeze(Ptild(k,:,:))*H'+R);
    
    %update state
    xhat(:,k) = xtild(:,k)+squeeze(K(k,:,:))*(z(:,k)-H*xtild(:,k));
    
    %update covariance
    P(k,:,:) = (I-squeeze(K(k,:,:))*H)*squeeze(Ptild(k,:,:));
   
end

%% Extract KF State Info
x_n_K = xhat(1,:);
x_e_K = xhat(2,:);
v_n_K = xhat(3,:);
v_e_K = xhat(4,:);

%% Calculations
sigma_xy_gps_emp = sqrt(std(x_m(:,1)-gps(:,1))^2+std(x_m(:,2)-gps(:,2))^2); %4.1
sigma_xy_KF_emp = sqrt(std(x_m(:,1)-x_n_K')^2+std(x_m(:,2)-x_e_K')^2); %4.2
sigma_xy_KF_P = sqrt(std(squeeze(P(:,1,1)))^2+std(squeeze(P(:,2,2)))^2);

innovation_all = zeros(size(gps'));
for i=1:n
    innovation_all(:,i) = gps(i,:)'-H*xtild(:,i); %4.3
    %innovation_stable = ;
end
%plots for innovation
figure(17);
plot(innovation_all(1,:));
hold on;
plot(innovation_all(2,:));

%% Calculate Errors
gps_err = n_gps;
kalman_err = xhat(1:2,:)-x_m';
v_kalman_err = xhat(3:4,:)-dx_m';

%% PLOTS REQUIRED: 2x1 subplot with position measurements (ref, meas, KF) and error (meas, KF)
%%Plotting Trajectories
figure(k);
%position
subplot(2,1,1);
plot(x_m(:,2),x_m(:,1));
hold on;
xlabel('East (m)');
ylabel('North (m)');
title('Trajectory Comparison of Raw GPS and Kalman Filter');
plot(gps(:,1),gps(:,2));
plot(x_e_K,x_n_K);
axis('equal');
legend('Ref','GPS','KF');

%%Plotting Errors
subplot(2,1,2);
scatter(gps_err(:,2),gps_err(:,1),'filled');
hold on;
scatter(kalman_err(2,:),kalman_err(1,:),'filled');
xbound = max(abs([kalman_err(1,:)';gps_err(:,1)]));
ybound = max(abs([kalman_err(2,:)';gps_err(:,2)]));
xlim([-xbound,xbound]);
ylim([-ybound,ybound]);
xlabel('East Error (m)');
ylabel('North Error (m)');
legend('GPS','KF');


%%               : 2x1 subplot with velocity measurements (ref, KF) and error (KF)
figure(k+1);
%velocity
subplot(2,1,1);
plot(dx_m(:,2),dx_m(:,1));
hold on;
xlabel('East (m/s)');
ylabel('North (m/s)');
title('Velocity Comparison of Kalman Filter and the Reference');
plot(v_e_K,v_n_K);
axis('equal');
legend('Ref','KF');

%%Plotting Errors
subplot(2,1,2);
scatter(v_kalman_err(2,:),v_kalman_err(1,:),'filled');
hold on;
xbound = max(abs([v_kalman_err(1,:)';gps_err(:,1)]));
ybound = max(abs([v_kalman_err(2,:)';gps_err(:,2)]));
xlim([-xbound,xbound]);
ylim([-ybound,ybound]);
xlabel('East Error (m/s)');
ylabel('North Error (m/s)');
legend('KF');

%%               : Histogram innovation histogram for north and east



