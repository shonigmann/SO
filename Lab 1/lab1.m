set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultLineLineWidth',1.5);

%%%% Lab 1
%% Part A
%1) Generate 3 random sequences (white noise) - calculations in Random Walk
%function
%2) Generate random walks using sequences
%sequences are denoted s, walks are denoted rw
[s1,rw1] = randomWalk(200000,1);
[s2,rw2] = randomWalk(200000,2);
[s3,rw3] = randomWalk(200000,3);

%for plotting
x = 1:200000;

%3) Generate 1st order guass-markov process for two correlation times
dt = 1;
%a) tau = 2000
t=2000;

fogm1_1 = firstOrderGaussMarkov(t,dt,s1);
fogm2_1 = firstOrderGaussMarkov(t,dt,s2);
fogm3_1 = firstOrderGaussMarkov(t,dt,s3);
%b) tau = 500
t=500;

fogm1_2 = firstOrderGaussMarkov(t,dt,s1);
fogm2_2 = firstOrderGaussMarkov(t,dt,s2);
fogm3_2 = firstOrderGaussMarkov(t,dt,s3);

S = [s1';s2';s3'];
RW = [rw1';rw2';rw3'];
FOGMa = [fogm1_1';fogm2_1';fogm3_1'];
FOGMb = [fogm1_2';fogm2_2';fogm3_2'];

%Save data in a text file
fileID = fopen('stochastic_process2000.txt','w');
fprintf(fileID,'%8.8f, %8.8f, %8.8f\n',[fogm1_1,fogm2_1,fogm3_1]);
fclose(fileID);

fileID = fopen('stochastic_process500.txt','w');
fprintf(fileID,'%8.8f, %8.8f, %8.8f\n',[fogm1_2,fogm2_2,fogm3_2]);
fclose(fileID);

fileID = fopen('rw.txt','w');
fprintf(fileID,'%8.8f, %8.8f, %8.8f\n',RW);
fclose(fileID);

fileID = fopen('wn.txt','w');
fprintf(fileID,'%8.8f, %8.8f, %8.8f\n',S);
fclose(fileID);

%sanity check:
figure();
title('Original Signals');
hold on;
subplot(4,1,1);
plot(x,s1,x,s2,x,s3);
legend('1','2','3');
ylabel('Random Sequence');
subplot(4,1,2);
plot(x,rw1,x,rw2,x,rw3);
legend('1','2','3');
ylabel('Random Walks');
subplot(4,1,3);
plot(x,fogm1_1,x,fogm2_1,x,fogm3_1);
ylabel('G-M P (T=2000)');
legend('1','2','3');
subplot(4,1,4);
plot(x,fogm1_2,x,fogm2_2,x,fogm3_2);
ylabel('G-M P (T=500)');
legend('1','2','3');
xlabel('sample number');

%% Part B
%4) compute the noise characteristics for each sequence by a:
%autocorrelation function b: Power spectral density, c: (optional) allan
%variance
for i=1:3
    
%a) autocorrelation function
    [C_S(:,i),L_S] = xcorr(S(i,:));
    [C_RW(:,i),L_RW] = xcorr(RW(i,:));
    [C_FOGMa(:,i),L_FOGMa] = xcorr(FOGMa(i,:));
    [C_FOGMb(:,i),L_FOGMb] = xcorr(FOGMb(i,:));
        
%b) power spectral density
    [H_S(:,i),w1] = pwelch(C_S(:,i));
    [H_RW(:,i),w2] = pwelch(C_RW(:,i));
    [H_FOGMa(:,i),w3] = pwelch(C_FOGMa(:,i));
    [H_FOGMb(:,i),w4] = pwelch(C_FOGMb(:,i));
    
%c) allan variance
    %%skipping for now... 
end

%plotting values from Part B:
figure(2);
x = 1:size(C_S,1);
for i=1:3
    subplot(2,1,1);
    hold on;
    plot(L_S,C_S(:,i));
end
title('White Noise');
legend('WN1','WN2','WN3');
ylabel('Autocorrelation Function');
xlabel('Lag Time');
ylabel('Autocorrelation Sequence');

figure(3);
for i=1:3
    subplot(2,1,1);
    hold on;
    plot(L_RW,C_RW(:,i));
end
title('Random Walk');
legend('RW1','RW2','RW3');
xlabel('Lag Time');
ylabel('Autocorrelation Sequence');

figure(4);
for i=1:3
    subplot(2,1,1);
    hold on;
    plot(L_FOGMa,C_FOGMa(:,i));
end
title('Gauss-Markov (T=2000)');
legend('FOGMa1','FOGMa2','FOGMa3');
xlabel('Lag Time');
ylabel('Autocorrelation Sequence');

figure(5);
for i=1:3
    subplot(2,1,1);
    hold on;
    plot(L_FOGMb,C_FOGMb(:,i));
end
title('Gauss-Markov (T=500)');
legend('FOGMb1','FOGMb2','FOGMb3');
xlabel('Lag Time');
ylabel('Autocorrelation Sequence');
   
x=w1;
figure(2);
for i=1:3
    subplot(2,1,2);
    hold on;
    plot(x,10*log10(H_S(:,i)));
end
title('White Noise');
legend('WN1','WN2','WN3');
ylabel('Power Spectral Density (dB)');
xlabel('Frequency (rad/sample)');
set(gca,'XScale','log');

figure(3);
for i=1:3
    subplot(2,1,2);
    hold on;
    plot(x,10*log10(H_RW(:,i)));
end
title('Random Walk');
legend('RW1','RW2','RW3');
ylabel('Power Spectral Density (dB)');
xlabel('Frequency (rad/sample)');
set(gca,'XScale','log');

figure(4);
for i=1:3
    subplot(2,1,2);
    hold on;
    plot(x,10*log10(H_FOGMa(:,i)));
end
title('Gauss-Markov (T=2000)');
legend('FOGMa1','FOGMa2','FOGMa3');
ylabel('Power Spectral Density (dB)');
xlabel('Frequency (rad/sample)');
set(gca,'XScale','log');

figure(5);
for i=1:3
    subplot(2,1,2);
    hold on;
    plot(x,10*log10(H_FOGMb(:,i)));
end
title('Gauss-Markov (T=500)');
legend('FOGMb1','FOGMb2','FOGMb3');
ylabel('Power Spectral Density (dB)');
xlabel('Frequency (rad/sample)');
set(gca,'XScale','log');

%5) verify graphically the changes of these characteristics for each type
%of noise. make a plot for each type of noise while coggling the colors 


%6) compare your findings with the values determined by the online tool for
%noise characterization via GMWM. Note: only need to upload one sequence
%per each noise



