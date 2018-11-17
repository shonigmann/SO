%% turning on MATLAB connector
connector on 123456 % MAKE SURE YOUR PHONE IS CONNECTED BEFORE PROCEEDING
%% creating mobiledev object and setting the sample rate
m = mobiledev; m.SampleRate = 100; disp(m)
%% enabling accelerometers and gyroscopes
m.AccelerationSensorEnabled = 1;
m.AngularVelocitySensorEnabled = 1;
disp(m)
%% logging the data, MAKE SURE YOUR PHONE IS AT REST
m.Logging = 1;
%pause(15);
pause(15*60 + 60); % time in seconds, with a margin for delays in logging
m.Logging = 0;
%% retreiving logged data
[a, t_a] = accellog(m); [w, t_w] = angvellog(m);
%% plots
figure(1); plot(t_a,a); xlabel('t [s]'); ylabel('accel. [m/s^2]')
figure(2); plot(t_w,w); xlabel('t [s]'); ylabel('rot. rate [rad/s]')