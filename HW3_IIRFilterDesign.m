%% Shifra Abittan
% Professor Fontaine
% DSP ECE310
% Problem Set 3: IIR Filter Design

%% Specifications

% Sampling Rate
f_samp = 400*(10^3);
w_samp = 2*pi*f_samp;

Nyquist = f_samp./2

% W (capital) = analog
% w (lowercase) = digital

% ANALOG: Convert frequency from Hz to rad/sec
W_pHI = 2*pi*130*(10^3);
W_pLO = 2*pi*100*(10^3);

W_sHI = 2*pi*140*(10^3);
W_sLO = 2*pi*90*(10^3);

% DIGITAL: Normalize frequency by dividing f by Nyquist bandwidth???

% DIDNT USE: Use the formula w = (2*pi*f)./f_samp to determine normalized
% digital frequency; units radian
%{
w_pHI = W_pHI./f_samp;
w_pLO = W_pLO./f_samp;

w_sHI = W_sHI./f_samp;
w_sLO = W_sLO./f_samp;
%}

w_pHI = 130*(10^3)./Nyquist;
w_pLO = 100*(10^3)./Nyquist;

w_sHI = 140*(10^3)./Nyquist;
w_sLO = 90*(10^3)./Nyquist;

% Attenuation Values
rs = 30;
rp = 2;

B = W_pHI - W_pLO;
w_0 = (W_pHI*W_pLO)^0.5;


%% Generate Analog Filters

% Determine which of the two frequencies to use
w_LPproto1_s = abs((W_sLO^2 - w_0^2)/(B*W_sLO));
w_LPproto2_s = abs((W_sHI^2 - w_0^2)/(B*W_sHI)); % use this for stopband freq
w_LPproto_p = 1;

% Analog Butterworth Order
[N_Butter_Anlg, Wn_b] = buttord([W_pLO,W_pHI],[W_sLO,W_sHI],rp,rs,'s');
% Analog Butterworth Filter
[Z_Butter_Anlg,P_Butter_Anlg,K_Butter_Anlg] = butter(N_Butter_Anlg,Wn_b,'s');
[B_butter_anlg,A_butter_anlg] = butter(N_Butter_Anlg,Wn_b,'s'); %for impinvar only

% Analog Chebychev1 Order
[N_Cheb1_Anlg, Wp] = cheb1ord([W_pLO,W_pHI],[W_sLO,W_sHI],rp,rs,'s');
% Analog Chebychev1 Filter
[Z_Cheb1_Anlg,P_Cheb1_Anlg,K_Cheb1_Anlg] = cheby1(N_Cheb1_Anlg,rp,Wp,'s');

% Analog Chebychev2 Order
[N_Cheb2_Anlg, Ws] = cheb2ord([W_pLO,W_pHI],[W_sLO,W_sHI],rp,rs,'s');
% Analog Chebychev2 Filter
[Z_Cheb2_Anlg,P_Cheb2_Anlg,K_Cheb2_Anlg] = cheby2(N_Cheb2_Anlg,rs,Ws,'s');

% Analog Elliptic Order
[N_Ellip_Anlg, Wp] = ellipord([W_pLO,W_pHI],[W_sLO,W_sHI],rp,rs,'s');
% Analog Elliptic Filter
[Z_Ellip_Anlg,P_Ellip_Anlg,K_Ellip_Anlg] = ellip(N_Ellip_Anlg,rp,rs,Wp,'s');

%%Multiply by 2 for all orders

%% Generate Digital Filters (via bilinear transform)

% Digital Butterworth Order
[N_Butter_Dig, Wn_b] = buttord([w_pLO w_pHI],[w_sLO w_sHI],rp,rs);
% Digital Butterworth Filter
[Z_Butter_Dig,P_Butter_Dig,K_Butter_Dig] = butter(N_Butter_Dig,Wn_b);

% Digital Chebychev1 Order
[N_Cheb1_Dig, Wp] = cheb1ord([w_pLO,w_pHI],[w_sLO,w_sHI],rp,rs);
% Digital Chebychev1 Filter
[Z_Cheb1_Dig,P_Cheb1_Dig,K_Cheb1_Dig] = cheby1(N_Cheb1_Dig,rp,Wp);

% Digital Chebychev2 Order
[N_Cheb2_Dig, Ws] = cheb2ord([w_pLO,w_pHI],[w_sLO,w_sHI],rp,rs);
% Digital Chebychev2 Filter
[Z_Cheb2_Dig,P_Cheb2_Dig,K_Cheb2_Dig] = cheby2(N_Cheb2_Dig,rs,Ws);

% Digital Elliptic Order
[N_Ellip_Dig, Wp] = ellipord([w_pLO,w_pHI],[w_sLO,w_sHI],rp,rs);
% Digital Elliptic Filter
[Z_Ellip_Dig,P_Ellip_Dig,K_Ellip_Dig] = ellip(N_Ellip_Dig,rp,rs,Wp);


%% Generate Digital Filters (via impulse invariance)

% Iminvar Filter from Analog Butterworth
% Find num and denom of analog transfer fn
[B_butt_anlg,A_butt_anlg] = zp2tf(Z_Butter_Anlg,P_Butter_Anlg,K_Butter_Anlg);
% Create filter
[B_butt_imp,A_butt_imp] = impinvar(B_butt_anlg,A_butt_anlg,f_samp);
% Order
N_Butter_Imp = filtord(B_butt_imp,A_butt_imp);
% Z,P,K form
[Z_Butter_Imp,P_Butter_Imp,K_Butter_Imp]= tf2zpk(B_butt_imp,A_butt_imp);

% Iminvar Filter from Analog Chebychev1
% Find num and denom of analog transfer fn
[B_Cheb1_anlg,A_Cheb1_anlg] = zp2tf(Z_Cheb1_Anlg,P_Cheb1_Anlg,K_Cheb1_Anlg);
% Create filter
[B_Cheb1_imp,A_Cheb1_imp] = impinvar(B_Cheb1_anlg,A_Cheb1_anlg,f_samp);
% Order
N_Cheb1_Imp = filtord(B_Cheb1_imp,A_Cheb1_imp);
% Z,P,K form
[Z_Cheb1_Imp,P_Cheb1_Imp,K_Cheb1_Imp]= tf2zpk(B_Cheb1_imp,A_Cheb1_imp);

% Iminvar Filter from Analog Chebychev2
% Find num and denom of analog transfer fn
[B_Cheb2_anlg,A_Cheb2_anlg] = zp2tf(Z_Cheb2_Anlg,P_Cheb2_Anlg,K_Cheb2_Anlg);
% Create filter
[B_Cheb2_imp,A_Cheb2_imp] = impinvar(B_Cheb2_anlg,A_Cheb2_anlg,f_samp);
% Order
N_Cheb2_Imp = filtord(B_Cheb2_imp,A_Cheb2_imp);
% Z,P,K form
[Z_Cheb2_Imp,P_Cheb2_Imp,K_Cheb2_Imp]= tf2zpk(B_Cheb2_imp,A_Cheb2_imp);

% Iminvar Filter from Analog Elliptic
% Find num and denom of analog transfer fn
[B_Ellip_anlg,A_Ellip_anlg] = zp2tf(Z_Ellip_Anlg,P_Ellip_Anlg,K_Ellip_Anlg);
% Create filter
[B_Ellip_imp,A_Ellip_imp] = impinvar(B_Ellip_anlg,A_Ellip_anlg,f_samp);
% Order
N_Ellip_Imp = filtord(B_Ellip_imp,A_Ellip_imp);
% Z,P,K form
[Z_Ellip_Imp,P_Ellip_Imp,K_Ellip_Imp]= tf2zpk(B_Ellip_imp,A_Ellip_imp);


%% Pole Zero Plots

% Analog
figure
sgtitle('Analog Filters Pole Zero Plots')

subplot(2,2,1)
zplane(Z_Butter_Anlg,P_Butter_Anlg)
title('Butterworth')

subplot(2,2,2)
zplane(Z_Cheb1_Anlg,P_Cheb1_Anlg)
title('Chebychev I')

subplot(2,2,3)
zplane(Z_Cheb2_Anlg,P_Cheb2_Anlg)
title('Chebychev II')

subplot(2,2,4)
zplane(Z_Ellip_Anlg,P_Ellip_Anlg)
title('Elliptic')


% Digital
figure
sgtitle('Digital Filters (bilinear transform) Pole Zero Plots')

subplot(2,2,1)
zplane(Z_Butter_Dig,P_Butter_Dig)
title('Butterworth')

subplot(2,2,2)
zplane(Z_Cheb1_Dig,P_Cheb1_Dig)
title('Chebychev I')

subplot(2,2,3)
zplane(Z_Cheb2_Dig,P_Cheb2_Dig)
title('Chebychev II')

subplot(2,2,4)
zplane(Z_Ellip_Dig,P_Ellip_Dig)
title('Elliptic')


% Digital via impulse invariance
figure
sgtitle('Digital Filters via Impulse Invariance Pole Zero Plots')

subplot(2,2,1)
zplane(Z_Butter_Imp,P_Butter_Imp)
title('ImpInVar from Butterworth')

subplot(2,2,2)
zplane(Z_Cheb1_Imp,P_Cheb1_Imp)
title('ImpInVar from Chebychev1')

subplot(2,2,3)
zplane(Z_Cheb2_Imp,P_Cheb2_Imp)
title('ImpInVar from Chebychev2')

subplot(2,2,4)
zplane(Z_Ellip_Imp,P_Ellip_Imp)
title('ImpInVar from Elliptic')

%% Magnitude Response Plots

% Frequency in Hz from DC to Nyquist bandwidth
F = 0:1:Nyquist;
% Convert frequency to W(rad/sec)
W = 2.*pi.*F;

figure
sgtitle('Analog Filters Magnitude Response')

%% Analog Butterworth
[b_Butter_Anlg,a_Butter_Anlg] = zp2tf(Z_Butter_Anlg,P_Butter_Anlg,K_Butter_Anlg);
[H_Butter_Anlg,W_B] = freqs(b_Butter_Anlg,a_Butter_Anlg,W);

subplot(2,2,1)
plot(F,20*log10(abs(H_Butter_Anlg)))

title('Analog Butterworth Filter: Magnitude Response')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

upperdb = 1;
lowerdb = -50;

% Filter specifications graphed
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


%% Analog Chebychev1
[b_Cheb1_Anlg,a_Cheb1_Anlg] = zp2tf(Z_Cheb1_Anlg,P_Cheb1_Anlg,K_Cheb1_Anlg);
[H_C1_Anlg,W_C1] = freqs(b_Cheb1_Anlg,a_Cheb1_Anlg,W);

subplot(2,2,2)
plot(F,mag2db(abs(H_C1_Anlg)))

title('Analog Chebychev1 Filter: Magnitude Response')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


%% Analog Chebychev2
[b_Cheb2_Anlg,a_Cheb2_Anlg] = zp2tf(Z_Cheb2_Anlg,P_Cheb2_Anlg,K_Cheb2_Anlg);
[H_C2_Anlg,W_C2] = freqs(b_Cheb2_Anlg,a_Cheb2_Anlg,W);

subplot(2,2,3)
plot(F,mag2db(abs(H_C2_Anlg)))

title('Analog Chebychev2 Filter: Magnitude Response')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


%% Analog Elliptic
[b_Ellip_Anlg,a_Elliptic_Anlg] = zp2tf(Z_Ellip_Anlg,P_Ellip_Anlg,K_Ellip_Anlg);
[H_Ellip_Anlg,W_Elliptic] = freqs(b_Ellip_Anlg,a_Elliptic_Anlg,W);

subplot(2,2,4)
plot(F,mag2db(abs(H_Ellip_Anlg)))

title('Analog Elliptic Filter: Magnitude Response')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


figure
sgtitle('Digital Filters Magnitude Response')

%% Digital Butterworth
[r,N] = size(F);

[b_Butter_Dig,a_Butter_Dig] = zp2tf(Z_Butter_Dig,P_Butter_Dig,K_Butter_Dig);
[H_Butter_Dig,F_Butter_Dig] = freqz(b_Butter_Dig,a_Butter_Dig,N,f_samp);

subplot(2,2,1)
plot(F,20*log10(abs(H_Butter_Dig)))

title('Butterworth')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

upperdb = 1;
lowerdb = -50;

% Filter specifications graphed
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')


%% Digital Chebychev1
[b_Cheb1_Dig,a_Cheb1_Dig] = zp2tf(Z_Cheb1_Dig,P_Cheb1_Dig,K_Cheb1_Dig);
[H_C1_Dig,F_C1_Dig] = freqz(b_Cheb1_Dig,a_Cheb1_Dig,N,f_samp);

subplot(2,2,2)
plot(F,20*log10(abs(H_C1_Dig)))

title('Chebychev1')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')

%% Digital Chebychev2
[b_Cheb2_Dig,a_Cheb2_Dig] = zp2tf(Z_Cheb2_Dig,P_Cheb2_Dig,K_Cheb2_Dig);
[H_C2_Dig,F_C2_Dig] = freqz(b_Cheb2_Dig,a_Cheb2_Dig,N,f_samp);

subplot(2,2,3)
plot(F,20*log10(abs(H_C2_Dig)))

title('Chebychev2')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


%% Digital Elliptic
[b_Ellip_Dig,a_Elliptic_Dig] = zp2tf(Z_Ellip_Dig,P_Ellip_Dig,K_Ellip_Dig);
[H_Ellip_Dig,F_Elliptic_Dig] = freqz(b_Ellip_Dig,a_Elliptic_Dig,N,f_samp);

subplot(2,2,4)
plot(F,20*log10(abs(H_Ellip_Dig)))

title('Elliptic')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


figure
sgtitle('Digital Filters via Impulse Invariance Magnitude Response')

%% Iminvar from Butterworth
[b_Butter_Imp,a_Butter_Imp] = zp2tf(Z_Butter_Imp,P_Butter_Imp,K_Butter_Imp);
[H_Butter_Imp,F_Butter_Imp] = freqz(b_Butter_Imp,a_Butter_Imp,N,f_samp);

subplot(2,2,1)
plot(F,20*log10(abs(H_Butter_Imp)))

title('Impulse Response Matching Butterworth Filter')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

upperdb = 1;
lowerdb = -50;

% Filter specifications graphed
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


%% Impinvar from Chebychev1
[b_Cheb1_Imp,a_Cheb1_Imp] = zp2tf(Z_Cheb1_Imp,P_Cheb1_Imp,K_Cheb1_Imp);
[H_C1_Imp,F_C1_Imp] = freqz(b_Cheb1_Imp,a_Cheb1_Imp,N,f_samp);

subplot(2,2,2)
plot(F,20*log10(abs(H_C1_Imp)))

title('Impulse Response Matching Chebychev1 Filter')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


%% Impinvar from Chebychev2
[b_Cheb2_Imp,a_Cheb2_Imp] = zp2tf(Z_Cheb2_Imp,P_Cheb2_Imp,K_Cheb2_Imp);
[H_C2_Imp,F_C2_Imp] = freqz(b_Cheb2_Imp,a_Cheb2_Imp,N,f_samp);

subplot(2,2,3)
plot(F,mag2db(abs(H_C2_Imp)))

title('Impulse Response Matching Elliptic Filter')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


%% Impinvar from Elliptic
[b_Ellip_Imp,a_Elliptic_Imp] = zp2tf(Z_Ellip_Imp,P_Ellip_Imp,K_Ellip_Imp);
[H_Ellip_Imp,F_Ellip_Imp] = freqz(b_Ellip_Imp,a_Elliptic_Imp,N,f_samp);

subplot(2,2,4)
plot(F,mag2db(abs(H_Ellip_Imp)))

title('Iminvar from Elliptic')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')
ylim([-50 1])
xlim([0 Nyquist])

% Plot specifications
hold on
plot(repmat(100*(10^3),1,(upperdb - lowerdb) + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','r')
hold on
plot(repmat(90*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(repmat(140*(10^3),1,upperdb - lowerdb + 1),(lowerdb:1:upperdb),'--','Color','g')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')
hold on
plot(F,repmat(-30,1,200001),':','Color','k')
%legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')


%% Order of Filters
% Multiply by 2 because these are bandpass filters.  The order given by
% MATLAB is for the underlying LP prototype.
Type_of_Filter = ['Analog' 'Digital' 'Digital via ImpInVar'].';
Butterworth_Order = [N_Butter_Anlg*2 N_Butter_Dig*2 N_Butter_Imp].'; %PROBLEM WITH 16 Butt DIg??
Chebychev1_Order = [N_Cheb1_Anlg*2 N_Cheb1_Dig*2 N_Cheb1_Imp].';
Chebychev2_Order = [N_Cheb2_Anlg*2 N_Cheb2_Dig*2 N_Cheb2_Imp].';
Elliptic_Order = [N_Ellip_Anlg*2 N_Ellip_Dig*2 N_Ellip_Imp].';

table(Butterworth_Order,Chebychev1_Order,Chebychev2_Order,Elliptic_Order,'RowNames',{'Analog' 'Digital' 'Digital via ImpInVar'})

%% Time Constants
% Analog Time Constants are determined by the pole closest to the jw axis.
% Because all of the poles are in the left hand plane (all of the filters
% are stable), to determine the time constant, identify the pole with the 
% largest real part. Then use this pole in the formula: 
% TimeConstant(tau) = mag(real part)

% Butterworth
tau_ButterAnalog = abs(max(real(P_Butter_Anlg)))

% Chebychev 1
tau_Cheb1Analog = abs(max(real(P_Cheb1_Anlg)))

% Chebychev 2
tau_Cheb2Analog = abs(max(real(P_Cheb2_Anlg)))

% Elliptic
tau_EllipticAnalog = abs(max(real(P_Ellip_Anlg)))


% Digital Time Constants are determined by the pole closest to the unit
% circle. Because all of the poles are inside the unit circle/stable
% filters, the pole corresponding to the time constant is the one with the
% largest magnitude. Then use this pole in the formula:
% TimeConstant(tau) = T/abs(ln(abs(p)))

T = 400*(10^3);
% Butterworth
P_BD = max(abs(P_Butter_Dig));
tau_ButterDig = T./abs(ln(abs(P_BD)))

% Chebychev 1
P_C1D = max(abs(P_Cheb1_Dig));
tau_Cheb1Dig = T./abs(ln(abs(P_C1D)))

% Chebychev 2
P_C2D = max(abs(P_Cheb2_Dig));
tau_Cheb2Dig = T./abs(ln(abs(P_C2D)))

% Elliptic
P_ED = max(abs(P_Ellip_Imp));
tau_EllipDig = T./abs(ln(abs(P_ED)))

% Iminvar derived from Butterworth
P_BI = max(abs(P_Butter_Imp));
tau_ButterIIV = T./abs(ln(abs(P_B1)))

% Iminvar derived from Chebychev 1
P_C1I = max(abs(P_Cheb1_Imp));
tau_Cheb1IIV = T./abs(ln(abs(P_C1I)))

% Iminvar derived from Chebychev 2
P_C2I = max(abs(P_Cheb2_Imp));
tau_Cheb2IIV = T./abs(ln(abs(P_C2I)))

% Iminvar derived from Elliptic
P_EI = max(abs(P_Ellip_Imp));
tau_EllipIIV = T./abs(ln(abs(P_E1)))
