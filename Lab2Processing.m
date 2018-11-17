%% Simon Honigmann
% Sensor Orientation Lab 2
% Group: 7
% Sensors: AIRINS/XSEA (Mechanical) & NavChip/Intersense (MEMS)
% Sensor Analyzed: Gyroscope
% Axes: Y & Z

%% Collect Data
clc;
load Sensor1and4Data; %saved relevant data as a mat file to not have to manually select every time

%extract and store sampling frequencies
fs_NC = 1/(dataNavChip(2,1)-dataNavChip(1,1));
fs_XS = 1/(dataXSEA(2,1)-dataXSEA(1,1));

%% Compute Characteristic Parameters

stdDevs = [std(dataNavChip(922:187222,3)), std(dataNavChip(922:187222,4));std(dataXSEA(:,3)),std(dataXSEA(:,4))];
stdDevs = [stdDevs,mean(stdDevs,2)];

means = [mean(dataNavChip(922:187222,3)), mean(dataNavChip(922:187222,4));mean(dataXSEA(:,3)),mean(dataXSEA(:,4))];
means = [means,mean(means,2)];

%% Plot Raw Data
if 0
figure(1);
subplot(2,1,1);
plot(dataNavChip(922:187222,1),dataNavChip(922:187222,3));
title('NavChip Gyroscope, Y Axis');
xlabel('time, (s)');
ylabel('Rotational Velocity (rad/s)');

subplot(2,1,2);
plot(dataNavChip(922:187222,1),dataNavChip(922:187222,4));
title('NavChip Gyroscope, Z Axis');
xlabel('time, (s)');
ylabel('Rotational Velocity (rad/s)');

figure(2);
subplot(2,1,1);
plot(dataXSEA(:,1),dataXSEA(:,3));
title('XSea Gyroscope, Y Axis');
xlabel('time, (s)');
ylabel('Rotational Velocity (rad/s)');

subplot(2,1,2);
plot(dataXSEA(:,1),dataXSEA(:,4));
title('XSea Gyroscope, Z Axis');
xlabel('time, (s)');
ylabel('Rotational Velocity (rad/s)');
end
%% Autocorrelation
[C_NCy,L_NCy] = xcorr(dataNavChip(922:187222,3),'coeff');
[C_NCz,L_NCz] = xcorr(dataNavChip(922:187222,4),'coeff');
[C_XSy,L_XSy] = xcorr(dataXSEA(:,3),'coeff');
[C_XSz,L_XSz] = xcorr(dataXSEA(:,4),'coeff');

figure(3);
subplot(2,1,1);
plot(L_NCy,(C_NCy));
hold on;
plot(L_NCz,(C_NCz));
title('Nav Chip Gyroscope Autocorrelation');
legend('Y axis', 'Z axis');
ylabel('Sample Autocorrelation Function');
xlabel('Lag Time (s)');

figure(4);
subplot(2,1,1);
plot(L_XSy,(C_XSy));
hold on;
plot(L_XSz,(C_XSz));
title('XSEA Gyroscope Autocorrelation');
legend('Y axis', 'Z axis');
ylabel('Sample Autocorrelation Function');
xlabel('Lag Time (s)');


% Power Spectral Density
[H_NCy,f_NCy] = pwelch(dataNavChip(922:187222,3),[],[],[],fs_NC);
[H_NCz,f_NCz] = pwelch(dataNavChip(922:187222,4),[],[],[],fs_NC);
[H_XSy,f_XSy] = pwelch(dataXSEA(:,3),[],[],[],fs_XS);
[H_XSz,f_XSz] = pwelch(dataXSEA(:,4),[],[],[],fs_XS);

figure(3);
subplot(2,1,2);
x=f_NCy;
y=10*log10(H_NCy);
x = [-fliplr(x')';x];
y = [fliplr(y')';y];
plot(x,y);
hold on;
x=f_NCz;
y=10*log10(H_NCz);
x = [-fliplr(x')';x];
y = [fliplr(y')';y];
plot(x,y);title('Nav Chip Gyroscope Power Spectral Density (dB/Hz)');
legend('Y axis', 'Z axis');
ylabel('Power Spectral Density');
xlabel('Frequency (Hz)');
set(gca,'XScale','log');

figure(4);
subplot(2,1,2);
x=f_XSy;
y=10*log10(H_XSy);
x = [-fliplr(x')';x];
y = [fliplr(y')';y];
plot(x,y);
hold on;
x=f_XSz;
y=10*log10(H_XSz);
x = [-fliplr(x')';x];
y = [fliplr(y')';y];
plot(x,y);
title('XSEA Gyroscope Power Spectral Density');
legend('Y axis', 'Z axis');
ylabel('Power Spectral Density (dB/Hz)');
xlabel('Frequency (Hz)');
set(gca,'XScale','log');

%% Allan Variance
if 0
    av_NCy = allandev(dataNavChip(922:187222,3),'NavChip Gyro Y Axis',5,'k',100);
    av_NCz = allandev(dataNavChip(922:187222,4),'NavChip Gyro Z Axis',6,'k',100);
    av_XSy = allandev(dataXSEA(:,3),'XSEA Gyro Y Axis',7,'k',200);
    av_XSz = allandev(dataXSEA(:,4),'XSEA Gyro Z Axis',8,'k',200);
end
%% Synthetic Noise for z axis of NavChip and Y axis of XSea
rng(1); %seed for repeatibility

%NavChip
std_z = stdDevs(1,2); %depends if I want to use std dev calculated for
%signal, or variance 
wn = randn(187222-922+1,4); %white noise process
[~,rw] = randomWalk(size(wn,1),12345);
%Quantized Noise process
f=100;
dt = 1/f;
Q = 4.52e-7; %assuming GMWM summary returns Q and not sigma^2
U=wn(:,1); %assuming the base noise process for quantized and GaussMarkov needs to be a different WN process for each
dU = (U(2:length(U))-U(1:(length(U)-1))); %approximation of derivative of noise, assuming that's what Udot is in equation 5.36
QN = sqrt(Q)*dU; %eq. from table 5.2

%Bias Instability
t = dataNavChip(922:187222,1)-dataNavChip(922,1);
Tbi = 100;
bi = biasInstability(t,wn(:,4),Tbi);

%Gauss-Markov Process
Beta1 = 1/687;
Beta2 = 1/(263e-3);
sigma2_1 = 1.38e-8;
sigma2_2 = 3.08e-7;

GM1 = (sigma2_1)*firstOrderGaussMarkov(1/Beta1,1/f,wn(1:186300,2));
GM2 = (sigma2_2)*firstOrderGaussMarkov(1/Beta2,1/f,wn(1:186300,3));

%Linear Combination of Processes Exhibited
NCNoiseZ = (bi(1:186300)+wn(1:186300,1))*std_z+QN;
%NCNoiseZ = QN+bi(1:186300)*std_z; %desperately trying to debug

%PLOTTING THINGS
%Autocorrelation
[C,L] = xcorr(NCNoiseZ,'coeff');
figure(3);
subplot(2,1,1);
plot(L,(C),'g');
hold on;
title('Simulated NavChip Gyroscope Autocorrelation');
ylabel('Sample Autocorrelation Function');
xlabel('Lag Time (s)');

% Power Spectral Density
[H,fp] = pwelch(NCNoiseZ,[],[],[],f);

subplot(2,1,2);
x=fp;
y=10*log10(H);
x = [-fliplr(x')';x];
y = [fliplr(y')';y];
plot(x,y,'g');
hold on;
title('Simulated NavChip Gyroscope PSD (dB/Hz)');
ylabel('Power Spectral Density');
xlabel('Frequency (Hz)');
set(gca,'XScale','log');

%allan variance
avNC_synth = allandev(NCNoiseZ,'NavChip Gyro Z Axis',1,'g',100);
%av_NCz = allandev(dataNavChip(922:187222,4),'NavChip Gyro Z Axis',9,'k',100);

%% XSea Gyro Synthetic:
if 0 %stop from graphing when not working on this section
    std_y = stdDevs(2,1);
    wn_y = std_y*randn(size(dataXSEA,1),1);

    %Autocorrelation
    [C_NCy,L_NCy] = xcorr(wn_y,'coeff');
    figure(4);
    subplot(2,1,1);
    plot(L_NCy,(C_NCy),'g');
    hold on;
    title('Simulated XSea Gyroscope Autocorrelation');
    ylabel('Sample Autocorrelation Function');
    xlabel('Lag Time (s)');

    % Power Spectral Density
    [H_NCy,f_NCy] = pwelch(wn_y,[],[],[],fs_NC);

    subplot(2,1,2);
    x=f_NCy;
    y=10*log10(H_NCy);
    x = [-fliplr(x')';x];
    y = [fliplr(y')';y];
    plot(x,y,'g');
    hold on;
    title('Simulated XSea Gyroscope PSD (dB/Hz)');
    ylabel('Power Spectral Density');
    xlabel('Frequency (Hz)');
    set(gca,'XScale','log');

    %allan variance
    av_XSy = allandev(dataXSEA(:,3),'XSea Gyro Y Axis',10,'k',200);
    av_XSy_synth = allandev(wn_y,'XSea Gyro Y Axis',10,'g',200);
end