%% Simon Honigmann
%  Sensor Orientation, Lab 6
%  Nov 9, 2018

%% Cleanup
clear
clc

graphs = 1; %1 to graph output, 0 to turn graphing off

%% Reference Values
mean_earth_rotation = 7.2921150e-5; %mean earths rotation, rad/s (NED)
grav_epfl = (980000+550)*10^-5; %gravity at epfl, m/s^2
%IMU Reported Values
lat_IMU = 46.52093816; %IMU reported latitude, deg
long_IMU = 6.56585746; %IMU reported longitude, deg
h_IMU = -200.983; %IMU reported heading, rad
r_IMU = -1.984*pi/180; %IMU reported roll, rad
p_IMU = 0.588*pi/180; %IMU reported pitch, rad (NED)
alt_IMU = 410.258; %IMU reported Altitude, m
speed_IMU = 0; %IMU reported speed norm, m/s

%% Lab Formatting: 
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultLineLineWidth',1.5);

%% Load and Preprocess Data for XSEA navigation sensor
load('lab6data.mat');
%data = readimu('1109_1407_PostProBinaryDecoded.imu');

t = data(:,1); %time
sample_rate = 1/(t(2)-t(1)); %Hz

start_time = 481565; %(s) taken from photo
[~, start_index] = min(abs(t-start_time));

%start_index = 489114; %found by looking at all data and seeing where my angle was... photo time didn't work for me... 
stop_index = floor(2*60*sample_rate+start_index);

t=t(start_index:stop_index);
gyro = data(start_index:stop_index,2:4); %x y and z measurements of gyro
accel = data(start_index:stop_index,5:7); %x y and z measurements of accelerometer

%% Averaging over the time interval:
gyro_avg = mean(gyro);
accel_avg = mean(accel);

%% Calculate the Norm of the Gyros
gyro_avg(:,4) = sqrt(gyro_avg(1).^2+gyro_avg(2).^2+gyro_avg(3).^2);
accel_avg(:,4) = sqrt(accel_avg(1).^2+accel_avg(2).^2+accel_avg(3).^2);

%% Plotting Raw Data
axisLabels = ['X','Y','Z','Norm'];
if graphs==1
    for i=1:3
        figure(1);
        subplot(3,1,i);
        plot(t,gyro(:,i));
        title(['XSea Gyroscope, ',axisLabels(i),' Axis']);
        xlabel('Time, (s)');
        ylabel('Angular Vel. (rad/s)');

        figure(2);
        subplot(3,1,i);
        plot(t,accel(:,i));
        title(['XSea Accelerometer, ',axisLabels(i),' Axis']);
        xlabel('Time, (s)');
        ylabel('Acceleration (m/s^2)');
    end
end
%% Converting Reference Frame
%IMU Uses Front Left Up = X Y Z
R = [1 0 0; 0 -1 0; 0 0 -1]; %transformation matrix from IMU to NED

gyro_avg(1:3) = R*gyro_avg(1:3)'; %convert from NWU to NED
accel_avg(1:3) = R*accel_avg(1:3)';%convert from NWU to NED

%% Calculate the Errors (Q2.1, Q3.1)
err_earth_rotation = gyro_avg(4)-mean_earth_rotation;
p_err_earth_rotation = (gyro_avg(4)-mean_earth_rotation)/mean_earth_rotation*100;
err_g = accel_avg(4)-grav_epfl;
p_err_g = (accel_avg(4)-grav_epfl)/grav_epfl*100;

tableQ21 = [gyro_avg(4),mean_earth_rotation,err_earth_rotation, p_err_earth_rotation]
tableQ34 = [accel_avg(4),grav_epfl,err_g, p_err_g]

%% Accelerometer Leveling
%Calculated Values
g = -accel_avg(4); 
r = asin(accel_avg(2)/g); % roll, rad
p = asin(-accel_avg(1)/g); % pitch, rad

%errors
err_r = r-r_IMU;
p_err_r = err_r/r_IMU*100;
err_p = p-p_IMU;
p_err_p = err_p/p_IMU*100;

tableQ4 = [r,r_IMU,err_r, p_err_r; ...
           p,p_IMU,err_p, p_err_p]

%creating local level frame (NED) from body frame (front left up):
Az = 0; %for now
R1 = [1     0       0; ...
      0 cos(r)   sin(r); ...
      0 -sin(r)  cos(r)];
  
R2 = [cos(p)  0 -sin(p); ...
      0       1     0;  ...
      sin(p)  0 cos(p)];
  
%rotating to local frame with arbitrary azimuth
omega_l = ((R1*R2)'*gyro_avg(1:3)')

%% Gyrocompassing
Az = atan2(omega_l(2),omega_l(1));
 if Az>0
     Az = Az-2*pi;
 end
Az_deg = Az*180/pi
Az_err = Az_deg - h_IMU

% properly determining the local frame now:
R3 = [cos(Az)  sin(Az) 0; ...
      -sin(Az)  cos(Az) 0; ...
      0         0      1];

%rotating to local frame with arbitrary azimuth
omega_l = (R3*R2*R1*gyro_avg(1:3)');

%% Determining Latitude
lat = acosd(omega_l(1)/gyro_avg(4)); %from eq: wN = wEcos(phi) on slide 1
err_lat = lat-lat_IMU;
p_err_lat = err_lat/lat_IMU*100; 
tableQ5 = [lat,lat_IMU,err_lat,p_err_lat]
