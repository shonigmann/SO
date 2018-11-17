%% Cleanup
clear
clc

%% Reference Values
mean_earth_rotation = 7.2921150e-5; %mean earths rotation, rad/s
grav_epfl = (980000+550)*10^-5; %gravity at epfl, m/s^2
%IMU Reported Values
lat_IMU = 46.52093816; %IMU reported latitude, deg
long_IMU = 6.56585746; %IMU reported longitude, deg
h_IMU = 200.983; %IMU reported heading, deg
r_IMU = -1.984*pi/180; %IMU reported roll, rad 
p_IMU = 0.588*pi/180; %IMU reported pitch, rad (CONVERTED TO NED)
alt_IMU = 410.258; %IMU reported Altitude, m
speed_IMU = 0; %IMU reported speed norm, m/s

%% Load and Preprocess Data for XSEA navigation sensor
load('lab6data.mat');
%data = readimu('1109_1407_PostProBinaryDecoded.imu');
num_points =10000;
step = floor((size(data,1)-10000)/num_points);
t_all = data(:,1); %time
sample_rate = 1/(t_all(2)-t_all(1)); %Hz


for index = 1:num_points

    %start_time = 481565; %(s) taken from photo
    %[~, start_index] = min(abs(t-start_time));

    %stop_index = floor(1.5*60*sample_rate+start_index);

    start_index = 1+(index-1)*step;
    stop_index = start_index + 10000; %start491203 end 501203

    if(index == 8470)
        tester = 0;
    end
    
    t=t_all(start_index:stop_index);
    gyro = data(start_index:stop_index,2:4); %x y and z measurements of gyro
    accel = data(start_index:stop_index,5:7); %x y and z measurements of accelerometer

    %% Averaging over the time interval:
    gyro_avg = mean(gyro);
    accel_avg = mean(accel);

    %% Calculate the Norm of the Gyros
    gyro_avg(:,4) = sqrt(gyro_avg(1).^2+gyro_avg(2).^2+gyro_avg(3).^2);
    accel_avg(:,4) = sqrt(accel_avg(1).^2+accel_avg(2).^2+accel_avg(3).^2);

    %% Converting Reference Frame
    %IMU Uses Front Left Up = X Y Z
    R = [1 0 0; 0 -1 0; 0 0 1]; %transformation matrix from IMU to NED

    gyro_avg(1:3) = R*gyro_avg(1:3)';
    accel_avg(1:3) = R*accel_avg(1:3)';

    %% Accelerometer Leveling
    %Calculated Values
    g = -accel_avg(4); %using equations that are not in the NED frame
    r = asin(accel_avg(2)/g); % roll, rad
    p = asin(accel_avg(1)/g); % pitch, rad

    %creating local level frame (NED) from body frame (front left up):
    Az = 0; %for now
    R1 = [1     0       0; ...
          0 cos(r)   -sin(r); ...
          0 sin(r)  cos(r)];

    R2 = [cos(p)  0 sin(p); ...
          0       1     0;  ...
          -sin(p)  0 cos(p)];

    R3 = [cos(Az)  -sin(Az) 0; ...
          sin(Az)  cos(Az) 0; ...
          0         0      1];

    %rotating to local frame with arbitrary azimuth
    omega_l = ((R1*R2)'*gyro_avg(1:3)');

    %% Gyrocompassing
    Az = atan2(omega_l(2),omega_l(1));
    Az_deg(index) = Az*180/pi;
    if(Az_deg(index)<0)
        Az_deg(index) = Az_deg(index) + 360;
    end

end

plot(Az_deg);