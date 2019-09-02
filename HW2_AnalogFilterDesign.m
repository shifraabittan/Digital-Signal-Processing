%% Shifra Abittan
% Professor Fontaine
% DSP ECE310
% Problem Set 2: Analog Filter Design

%% (A)

% Use MATLAB as a calculator

w_pHI = 2*pi*130*(10^3);
w_pLO = 2*pi*100*(10^3);

w_sHI = 2*pi*140*(10^3);
w_sLO = 2*pi*90*(10^3);

B = w_pHI - w_pLO
% B also equals 6*pi*(10^4)

w_0 = (w_pHI*w_pLO)^0.5
% w_0 also equals 2*pi*(10^5)*(1.3)^0.5

w_LPproto1 = abs((w_sLO^2 - w_0^2)/(B*w_sLO))
w_LPproto2 = abs((w_sHI^2 - w_0^2)/(B*w_sHI))

%% (B) Filter Order

% Use MATLAB as a calculator

rs = 30;
rp = 2;
ws = w_LPproto2;
wp = 1;

n_Butterworth = (0.5*log10(((10^(rs/10))-1)/((10^(rp/10))-1)))/log10(ws/wp)
n_Chebyshev = acosh((((10^(rs/10))-1)/((10^(rp/10))-1))^0.5)/acosh(ws/wp)

%% (C) Build filters directly

% Butterworth Order
[N_Butterworth, Wn_b] = buttord([w_pLO,w_pHI],[w_sLO,w_sHI],rp,rs,'s');
% Butterworth Filter
[Z_Butterworth,P_Butterworth,K_Butterworth] = butter(N_Butterworth,Wn_b,'s');
% Check butterworth filter order matches
% {function that calculates order always calculates underlying order (even
% though bandpass doubles the n, the buttord function returns the n of the
% low pass prototype)}
isequal(9,N_Butterworth)

% Chebychev1 Order
[N_Cheb1, Wp] = cheb1ord([w_pLO,w_pHI],[w_sLO,w_sHI],rp,rs,'s');
% Chebychev1 Filter
[Z_Cheb1,P_Cheb1,K_Cheb1] = cheby1(N_Cheb1,rp,Wp,'s');
% Check chebychev1 filter order matches
isequal(5,N_Cheb1)

% Chebychev2 Order
[N_Cheb2, Ws] = cheb2ord([w_pLO,w_pHI],[w_sLO,w_sHI],rp,rs,'s');
% Chebychev2 Filter
[Z_Cheb2,P_Cheb2,K_Cheb2] = cheby2(N_Cheb2,rs,Ws,'s')
% Check chebychev1 filter order matches
isequal(5,N_Cheb2)

% Elliptic Order
[N_Elliptic, Wp] = ellipord([w_pLO,w_pHI],[w_sLO,w_sHI],rp,rs,'s');
% Elliptic Filter
[Z_Elliptic,P_Elliptic,K_Elliptic] = ellip(N_Elliptic,rp,rs,Wp,'s');

%% (D) Zplane Filter Plots
figure
subplot(2,2,1)
zplane(Z_Butterworth,P_Butterworth)
title('Butterworth Filter: Pole Zero Plot')

subplot(2,2,2)
zplane(Z_Cheb1,P_Cheb1)
title('Chebychev1 Filter: Pole Zero Plot')

subplot(2,2,3)
zplane(Z_Cheb2,P_Cheb2)
title('Chebychev2 Filter: Pole Zero Plot')

subplot(2,2,4)
zplane(Z_Elliptic,P_Elliptic)
title('Elliptic Filter: Pole Zero Plot')

%% (E) Elliptic Lowpass Prototype
% Elliptic Order of LP prototype filter
[N_Elliptic_LPP, Wp_LPP] = ellipord(wp,ws,rp,rs,'s');
% Elliptic Filter of LP prototype filter
[Z_Elliptic_LPP,P_Elliptic_LPP,K_Elliptic_LPP] = ellip(N_Elliptic_LPP,rp,rs,Wp_LPP,'s')

figure
zplane(Z_Elliptic_LPP,P_Elliptic_LPP)
title('Low Pass Prototype Elliptic Filter: Pole Zero Plot')

%% (F) Magnitude and Phase Response Plots

% Frequency from DC to 200kHz
F = 0:1:200*(10^3);
% Convert frequency to W(rad/sec)
W = 2.*pi.*F;

%% Butterworth
[b_Butterworth,a_Butterworth] = zp2tf(Z_Butterworth,P_Butterworth,K_Butterworth);
[H_B,W_B] = freqs(b_Butterworth,a_Butterworth,W);

figure
% Butterworth Magnitude Plot
subplot(2,1,1)
plot(F,mag2db(abs(H_B)))
title('Butterworth Filter: Magnitude Response')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')

% Vertical Scale
lowerdb = -60;
upperdb = 0;

ylim([lowerdb upperdb])
xlim([0 200*(10^3)])

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
legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')

% Butterworth Frequency Plot
subplot(2,1,2)
plot(F,unwrap(angle(H_B))*180/pi)
title('Butterworth Filter: Phase Response')
xlabel('Frequency (Hz)') 
ylabel('Unwrapped Phase Response (degrees)')
xlim([0 200*(10^3)])

%% Chebychev1
[b_Cheb1,a_Cheb1] = zp2tf(Z_Cheb1,P_Cheb1,K_Cheb1);
[H_C1,W_C1] = freqs(b_Cheb1,a_Cheb1,W);

figure
% Chebychev1 Magnitude Plot
subplot(2,1,1)
plot(F,mag2db(abs(H_C1)))
title('Chebychev1 Filter: Magnitude Response')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')

ylim([lowerdb upperdb])
xlim([0 200*(10^3)])

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
legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')

% Chebychev1 Frequency Plot
subplot(2,1,2)
plot(F,unwrap(angle(H_C1))*180/pi)
title('Chebychev1 Filter: Phase Response')
xlabel('Frequency (Hz)') 
ylabel('Unwrapped Phase Response (degrees)')
xlim([0 200*(10^3)])

%% Chebychev1
[b_Cheb2,a_Cheb2] = zp2tf(Z_Cheb2,P_Cheb2,K_Cheb2);
[H_C2,W_C2] = freqs(b_Cheb2,a_Cheb2,W);

figure
% Chebychev2 Magnitude Plot
subplot(2,1,1)
plot(F,mag2db(abs(H_C2)))
title('Chebychev2 Filter: Magnitude Response')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')

ylim([lowerdb upperdb]);
xlim([0 200*(10^3)]);

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
legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')

% Chebychev2 Frequency Plot
subplot(2,1,2)
plot(F,unwrap(angle(H_C2))*180/pi)
title('Chebychev2 Filter: Phase Response')
xlabel('Frequency (Hz)') 
ylabel('Unwrapped Phase Response (degrees)')
xlim([0 200*(10^3)])

%% Elliptic
[b_Elliptic,a_Elliptic] = zp2tf(Z_Elliptic,P_Elliptic,K_Elliptic);
[H_Elliptic,W_Elliptic] = freqs(b_Elliptic,a_Elliptic,W);

figure
% Elliptic Magnitude Plot
subplot(2,1,1)
plot(F,mag2db(abs(H_Elliptic)))
title('Elliptic Filter: Magnitude Response')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')

ylim([lowerdb upperdb])
xlim([0 200*(10^3)])

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
legend('Magnitude Response','passband edge','passband edge','stopband edge','stopband edge','rp','rs')

% Elliptic Frequency Plot
subplot(2,1,2)
plot(F,unwrap(angle(H_Elliptic))*180/pi)
title('Elliptic Filter: Phase Response')
xlabel('Frequency (Hz)') 
ylabel('Unwrapped Phase Response (degrees)')
xlim([0 200*(10^3)])

%%MAYBE SCALE EACH VERTICAL AXIS SEPERATLY

%% (G) Zoomed in Magnitude Response

figure

% Plot Specs
plot(repmat(100*(10^3),1,(5 - lowerdb) + 1),(lowerdb:1:5),'--','Color','r')
hold on
plot(repmat(130*(10^3),1,5 - lowerdb + 1),(lowerdb:1:5),'--','Color','r')
hold on
plot(F,repmat(0,1,200001),':','Color','k')
hold on
plot(F,repmat(-2,1,200001),':','Color','k')

% Label and resize graph
title('Magnitude Response in the Passband for Various Filters')
xlabel('Frequency (Hz)') 
ylabel('Magnitude Response (dB)')

xlim([95*(10^3) 135*(10^3)])
ylim([-2.5 0.5])

% Butterworth Magnitude Plot
hold on
plot(F,mag2db(abs(H_B)))

% Chebychev1 Magnitude Plot
hold on
plot(F,mag2db(abs(H_C1)))

% Chebychev2 Magnitude Plot
hold on
plot(F,mag2db(abs(H_C2)))

% Elliptic Magnitude Plot
hold on
plot(F,mag2db(abs(H_Elliptic)))


legend('Passband Edge','Passband Edge','DC','rp','Butterworth Filter','Chebychev1 Filter','Chebychev2 Filter','Elliptic Filter')

%% (H)
B_db = mag2db(abs(H_B));
Butterworth_sedge1 = B_db(90*(10^3) + 1)
Butterworth_sedge2 = B_db(140*(10^3) + 1)

C1_db = mag2db(abs(H_C1));
Cheb1_sedge1 = C1_db(90*(10^3) + 1)
Cheb1_sedge2 = C1_db(140*(10^3) + 1)

C2_db = mag2db(abs(H_C2));
Cheb2_sedge1 = C2_db(90*(10^3) + 1)
Cheb2_sedge2 = C2_db(140*(10^3) + 1)

Elliptic_db = mag2db(abs(H_Elliptic));
Elliptic_sedge1 = Elliptic_db(90*(10^3) + 1)
Elliptic_sedge2 = Elliptic_db(140*(10^3) + 1)