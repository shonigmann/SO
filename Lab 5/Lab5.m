%% Simon Honigmann
% Sensor Orientation
% Lab 5: Inertial Navigation in 2D / Realistic Signal
% 11/02/18

%% Cleanup
clc;
clear;
close all;

%% Lab Formatting: 
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultLineLineWidth',1.5);

%% Definitions:
% m = 2D mapping frame (non-accelerating, non-rotating)
%       x_m1 = north
%       x_m2 = east
% m frame polar coordinates
%       psi = polar angle
%       r = radius

% b = body frame
        % alpha = azimuth (yaw) angle
        % x_b1 = tangential direction
        % x_b2 = radial direction (inwards)

%% Constants
r = 500; % Circle radius: 500 m, Angular speed ? = ?/100;
omega = pi/100; % angular speed = pi/100
g = 9.81; %gravity, m/s^2

%only for 100Hz and Trapezoidal
sampling_rate = 100;
method = 'Trapezoidal';

f = sampling_rate;
sample_rate_text = num2str(sampling_rate,3);
num_rotations = 1; %number of time the vehicle goes  around the circle. increase to better show divergence

%% Assumptions/Constraints
dr = 0; %constant radius of motion in m frame
ddr = 0;
dpsi = omega; %constant angular velocity in m frame
ddpsi = 0;

%% Initial Conditions
psi_0 = 0; % Initial position: on North axis
alpha_0 = psi_0 + pi/2; 
x_0 = [r,0];
dx_0 = [0,omega*r]; % Initial velocity: north-axis: 0, east-axis: ? · radius

required_time = 2*num_rotations*pi/omega; %time for 1 full rotation [s]
t = (0:1/sampling_rate:required_time)'; %time vector, [s]
samples = (1:length(t))';

%% Noise Definitions
b_g_spec = 10; %gyro bias non-SI, deg/hr
b_g = b_g_spec*(pi/180)/(3600)*ones(length(t),1); %SI gyro bias, rad/s

sigma_g_gm_spec = 0.005; %spec'd gyro correlated noise (1st order GM); deg/s/sqrt(Hz)
sigma_g_gm_SI = sigma_g_gm_spec*pi/180; %gyro correlated noise std. dev, rad/s/sqrt(Hz)

T_g = 100; %gyro correlation period, s
beta_g = 1/T_g; %gyro correlation frequency, Hz

sigma_g_wn_spec = 0.1; %gyro white noise strenght, deg/sqrt(hr)
sigma_g_wn_SI = sigma_g_wn_spec*(pi*sqrt(f)/(180*60));

b_a_spec = 1/1000; % accelerometer bias, g
b_a = b_a_spec*g*ones(length(t),1); %accelerometer bias, m/s^2

sigma_a_wn_spec = 50e-6; %accelerometer white noise strength, g/sqrt(Hz)
sigma_a_wn_SI = sigma_a_wn_spec*g*sqrt(f); %accelerometer white noise strength, m/s^2/sqrt(Hz)

var_checks = [b_g(1);sigma_g_gm_SI;T_g;sigma_g_wn_SI;b_a(1);sigma_a_wn_SI];


%Scale PSD values of sigma -> (nvm... done directly in noise generation)
 sigma_g_gm = sigma_g_gm_SI;
 sigma_g_wn = sigma_g_wn_SI;
 sigma_a_wn = sigma_a_wn_SI;

%% Generating Stochastic Process Components

seed = 1234567;
n = length(t);

[n_g_wn,~] = whiteNoiseRandomWalk(n,seed); %generate white noise process
n_g_wn = n_g_wn*sigma_g_wn;%gyro white noise

[n_a_wn,~] = whiteNoiseRandomWalk(n,seed+1); %accelerometer white noise
n_a_wn = n_a_wn*sigma_a_wn;

[n_wn,~] = whiteNoiseRandomWalk(n,seed+2); %one final random white noise process

n_wn = n_wn*sqrt(sigma_g_gm^2*(1-exp(-2*beta_g*1/sampling_rate)));% eq 32

%this is probably wrong... would be worth double checking... 
n_g_gm = firstOrderGaussMarkov(T_g,1/sampling_rate,n_wn); %1st order gauss markov process 
[H,f]=pwelch(n_g_gm,[],[],[],sampling_rate);

%selected combinations of noise components to include:
errors = [0,0,0,0,0; ... no error
          1,0,0,0,0; ... gyro bias only
          0,1,0,0,0; ... gyro gauss markov only
          0,0,1,0,0; ... gyro white noise only
          0,1,1,0,0; ... gyro white noise and gauss markov
          1,1,1,0,0; ... all gyro errors
          0,0,0,1,0; ... accel bias only
          0,0,0,0,1; ... accel white noise only
          0,0,0,1,1; ... all accel errors
          1,1,1,1,1]; %  all errors combined
for k=10:10
    %loop through all desired iterations:
    has_g_bias = errors(k,1);
    has_a_bias = errors(k,4);
    has_g_wn = errors(k,3);
    has_a_wn = errors(k,5);
    has_g_gm = errors(k,2);

    % Net Noise
    n_g = has_g_bias*b_g + has_g_wn*n_g_wn + has_g_gm*n_g_gm; %combined error for gyro
    n_a = has_a_bias*b_a + has_a_wn*n_a_wn; % combined error for accelerometer

    % Plotting as sanity check
    figure(k+10);
    subplot(2,2,1);
    plot(t,n_g_wn*has_g_wn);
    hold on;
    plot(t,n_g_gm*has_g_gm);
    plot(t,b_g*has_g_bias);
    
    title('Gyro Measurement Noise');
    xlabel('Time [s]');
    ylabel('Gyro Error [rad/s]');
    legend('White Noise','Gauss-Markov','Bias');

    figure(k+10);
    subplot(2,2,2);
    plot(t,n_g);

    title('Net Gyro Measurement Error');
    xlabel('Time [s]');
    ylabel('Gyro Error [rad/s]');

    subplot(2,2,3);
    plot(t,n_a_wn*has_a_wn);
    hold on;
    plot(t,b_a*has_a_bias);
    title({'Accelerometer Measurement Noise',''});
    xlabel('Time [s]');
    ylabel('Accelerometer Error [m/s^2]');
    legend('White Noise','Bias');

    subplot(2,2,4);
    plot(t,n_a);
    title({'Net Accelerometer Measurement Error',''});
    xlabel('Time [s]');
    ylabel('Accelerometer Error [m/s^2]');

    %% Generating Actual Trajectories
    psi = (psi_0:num_rotations*2*pi/(length(t)-1):num_rotations*2*pi)'; %polar angle, [rad]
    alpha = psi + pi/2; %azimuth angle, Initial azimuth: 90? (towards x-accelerometer)

    %use this line instead to have azimuth angle wrap at 2pi
    %%alpha = mod(psi + pi/2,2*pi); %azimuth angle, Initial azimuth: 90? (towards x-accelerometer)

    x_m = [r*cos(psi),r*sin(psi)]; % map frame, [north, east]
    dx_m = [dr.*cos(psi)-r.*sin(psi).*dpsi , r.*cos(psi).*dpsi + dr.*sin(psi)];

    %(definitely could have simplified this...)
    ddx_m = [ddr.*cos(psi) - dr.*sin(psi).*dpsi - (dr.*sin(psi).*dpsi + r.*cos(psi).*dpsi.^2 + r.*sin(psi).*ddpsi), ... %ddx_m1         
             dr.*cos(psi).*dpsi - r.*sin(psi).*dpsi.^2 + r.*cos(psi).*ddpsi + ddr.*sin(psi)+dr.*cos(psi).*dpsi ]; %ddx_m2

    dalpha = dpsi*ones(length(t),1); %deriving alpha wrt time = deriving psi wrt time 

    v_b = [ones(length(t),1)*omega*r,zeros(length(t),1)]; % got lazy... this works here because omega and r are constant. Would need to change this line if that was not the case
    a_b = [zeros(length(t),1),ones(length(t),1)*omega^2*r]; %got lazy... this works here because omega and r are constant. Would need to change this line if that was not the case

    %% Adding Noise to Sensors

    gyro = dalpha + n_g; %ideal gyroscope
    accel = a_b + n_a; %ideal accelerometer, f

    % Considering the initial conditions to be known, apply strapdown inertial navigation to
    % calculate the trajectory parameters (i.e. the attitude (azimuth) and 2D velocity and
    % position vectors, respectively).
    %[Rmb,Rbm] = RotationMatrix(alpha); %define rotation matrices for all angles, alpha

    %initialize variables with initial conditions
    alpha_sd = zeros(length(alpha),1); %strapdown azimuth
    alpha_sd(1) = alpha_0;
    v_sd = zeros(length(t),2); %strapdown velocity, map frame
    v_sd(1,:) = dx_0;
    x_sd = v_sd;
    x_sd(1,:) = x_0;

    %Localizing along trajectory. doing this with a loop. we'll see if i regret it later
    for i=2:length(t)
        if strcmp(method,'Rectangular')
            alpha_sd(i) = alpha_sd(i-1)+gyro(i)*(t(i)-t(i-1)); %strapdown azimuth
            [~,Rbm] = RotationMatrix(alpha_sd(i));
            Rbm = squeeze(Rbm);
            v_sd(i,:) = (v_sd(i-1,:)'+Rbm*accel(i,:)'*(t(i)-t(i-1)))';
            x_sd(i,:) = (x_sd(i-1,:)'+v_sd(i,:)'*(t(i)-t(i-1)))';
        elseif strcmp(method, 'Trapezoidal')
            alpha_sd(i) = alpha_sd(i-1)+1/2*(gyro(i)+gyro(i-1))*(t(i)-t(i-1));
            [~,Rbm] = RotationMatrix(alpha_sd(i)); %current Rotation Matrix
            Rbm = squeeze(Rbm);
            [~,Rbm_p] = RotationMatrix(alpha_sd(i-1)); %prior Rotation Matrix
            Rbm_p = squeeze(Rbm_p);
            v_sd(i,:) = (v_sd(i-1,:)'+1/2*(Rbm*accel(i,:)'+Rbm*accel(i-1,:)')*(t(i)-t(i-1)))';
            x_sd(i,:) = (x_sd(i-1,:)'+1/2*(v_sd(i,:)+v_sd(i-1,:))'*(t(i)-t(i-1)))';
        else 
            break %err. method invalid
        end
    end

    %% Plotting Trajectories
    figure(k);
    %position
    subplot(2,2,1);
    plot(x_m(:,2),x_m(:,1));
    hold on;
    xlabel('X_m_2, East (m)');
    ylabel('X_m_1, North (m)');
    title(['True vs Estimate Position (',method,' Integration, ',sample_rate_text,'Hz)']);
    plot(x_sd(:,2),x_sd(:,1));
    axis('equal');
    legend('True',strcat(method,', ',sample_rate_text,'Hz'));

    %% Calculate Errors
    err_x = x_sd-x_m;
    err_alpha = alpha_sd - alpha;
    err_v = v_sd - dx_m;

    err_x = [err_x,sqrt(err_x(:,1).^2+err_x(:,2).^2)]; %add third column for position error magnitude
    err_v = [err_v,sqrt(err_v(:,1).^2+err_v(:,2).^2)]; %add third column for position error magnitude

    err_table(k,:) = [mean(abs(err_alpha)),max(abs(err_alpha)),mean(abs(err_x(:,3))),max(abs(err_x(:,3))),mean(abs(err_v(:,3))),max(abs(err_v(:,3)))];

    %% Plotting Errors
    %position
    subplot(2,2,3);
    plot(t,err_x(:,1));
    hold on;
    plot(t,err_x(:,2));
    plot(t,err_x(:,3));
    xlabel('Time (s)');
    ylabel('X_m Error (m)');
    title(['Position Error Over Time (',method,' Integration, ',sample_rate_text,'Hz)']);
    legend('North','East','Error Magnitude');

    %azimuth
    subplot(2,2,2);
    plot(t,err_alpha);
    hold on;
    xlabel('Time (s)');
    ylabel('Azimuth Error (rad)');
    title({['Azimuth Angle Error Over Time (',method,' Integration, ',sample_rate_text,'Hz)'],' '});

    %velocity v1_m
    subplot(2,2,4);
    plot(t,err_v(:,1));
    hold on;
    plot(t,err_v(:,2));
    plot(t,err_v(:,3));
    xlabel('Time (s)');
    ylabel('V_m Error (m/s)');
    if(sampling_rate > 10)        
        title({['Velocity Error Over Time (',method,' Integration, ',sample_rate_text,'Hz)'],' '});
    else
        title({['Velocity Error Over Time (',method,' Integration, ',sample_rate_text,'Hz)']});
    end
    legend('North','East','Error Magnitude');
end
    
%% Functions
function [Rmb,Rbm] = RotationMatrix(alpha)

    Rmb = zeros(length(alpha),2,2); %preallocate memory
    Rbm = Rmb; 

    for i=1:length(alpha)
        Rmb(i,:,:) = [cos(alpha(i)),sin(alpha(i)); -sin(alpha(i)),cos(alpha(i))];%transformation matrix R(m->b)
        Rbm(i,:,:) = squeeze(Rmb(i,:,:))'; %transformation matrix R(b->m)
    end
    
end