%% Shifra Abittan
% Professor Fontaine
% DSP ECE310
% Problem Set 8: Quantization

%% 4. 
% Attenuation Values
rs = 30;
rp = 1.5;

% 8th Order Bandpass Elliptic Digital Filter
[Z,P,K] = ellip(4,rp,rs,[0.3 0.6]);
% The order of the filter is 8 but n in the ellip function refer to the
% order of the underlying analog lowpass prototype, which is (n_bandpass)/2

%% (a) Graph the magnitude response of transfer function 
[b,a] = zp2tf(Z,P,K);
[H,W] = freqz(b,a,500);

figure
plot(W./pi,20*log10(abs(H))) %confirm this freq normalization makes sense
title('Magnitude Response of 8th Order Bandpass Digital Elliptic Filter')
xlabel('Frequency') 
ylabel('Magnitude Response (dB)')

%% (b) SOS form with L-inf scaling 
[sos_up,g_up] = zp2sos(Z,P,K,'up','inf');
[sos_dwn,g_dwn] = zp2sos(Z,P,K,'down','inf');

%% (c)
% SOS UP
%{
% Compute transfer function of each second order stage
stage1_up = tf(sos_up(1,1:3),sos_up(1,4:6),-1);
stage2_up = tf(sos_up(2,1:3),sos_up(2,4:6),-1);
stage3_up = tf(sos_up(3,1:3),sos_up(3,4:6),-1);
stage4_up = tf(sos_up(4,1:3),sos_up(4,4:6),-1);
%}

% Compute poles and zeros of each stage
[z1_up,p1_up,k1_up] = tf2zpk(sos_up(1,1:3),sos_up(1,4:6));
[z2_up,p2_up,k2_up] = tf2zpk(sos_up(2,1:3),sos_up(2,4:6));
[z3_up,p3_up,k3_up] = tf2zpk(sos_up(3,1:3),sos_up(3,4:6));
[z4_up,p4_up,k4_up] = tf2zpk(sos_up(4,1:3),sos_up(4,4:6));

% Check that magnitude of poles and zeros are increasing:
% Each zero and pole magnitude will print out twice because they all appear
% as complex conjugates
z1_up_mag = abs(z1_up);
z2_up_mag = abs(z2_up);
z3_up_mag = abs(z3_up);
z4_up_mag = abs(z4_up);
%Zeros all have same magnitude of 1

p1_up_mag = abs(p1_up);
p2_up_mag = abs(p2_up);
p3_up_mag = abs(p3_up);
p4_up_mag = abs(p4_up);
%Poles increase in magnitude as advance from stage to stage

if (z1_up_mag <= z2_up_mag <= z3_up_mag <= z4_up_mag)
    disp('Zeros increase in magnitude (or remain the same) as advance from stage to stage for SOS-up, as expected')
end
%The above statement doesnt print out because all of the zeros are 1 and
%there is some rounding error. This rounding error means that the zeros are
%not exactly equal in magnitude so the if statement fails and returns 0.

if (p1_up_mag <= p2_up_mag <= p3_up_mag <= p4_up_mag)
    disp('Poles increase in magnitude as advance from stage to stage for SOS-up, as expected')
end

% SOS DOWN
% Compute transfer function of each second order stage
%{
stage1_dwn = tf(sos_dwn(1,1:3),sos_dwn(1,4:6),-1); 
stage2_dwn = tf(sos_dwn(2,1:3),sos_dwn(2,4:6),-1);
stage3_dwn = tf(sos_dwn(3,1:3),sos_dwn(3,4:6),-1);
stage4_dwn = tf(sos_dwn(4,1:3),sos_dwn(4,4:6),-1);
%}

% Compute poles and zeros of each stage
[z1_dwn,p1_dwn,k1_dwn] = tf2zpk(sos_dwn(1,1:3),sos_dwn(1,4:6));
[z2_dwn,p2_dwn,k2_dwn] = tf2zpk(sos_dwn(2,1:3),sos_dwn(2,4:6));
[z3_dwn,p3_dwn,k3_dwn] = tf2zpk(sos_dwn(3,1:3),sos_dwn(3,4:6));
[z4_dwn,p4_dwn,k4_dwn] = tf2zpk(sos_dwn(4,1:3),sos_dwn(4,4:6));

if (abs(z1_dwn) <= abs(z2_dwn) <= abs(z3_dwn) <= abs(z4_dwn))
    disp('Zeros decrease in magnitude (or remain the same) as advance from stage to stage for SOS-down, as expected')
end

if (abs(p1_up) <= abs(p2_up) <= abs(p3_up) <= abs(p4_up))
    disp('Poles decrease in magnitude as advance from stage to stage for SOS-down, as expected')
end

%% (d)
% UP
B1u = sos_up(1,1:3);
B2u = sos_up(2,1:3);
B3u = sos_up(3,1:3);
B4u = sos_up(4,1:3);
A1u = sos_up(1,4:6);
A2u = sos_up(2,4:6);
A3u = sos_up(3,4:6);
A4u = sos_up(4,4:6);

[H1u,W] = freqz(g_up,A1u,500); %g / A1
[H2u,W] = freqz(g_up*B1u,conv(A1u,A2u),500); %g*B1 / A1*A2
[H3u,W] = freqz(g_up*conv(B1u,B2u),conv(conv(A1u,A2u),A3u),500); %g*B1*B2 / A1*A2*A3
[H4u,W] = freqz(g_up*conv(conv(B1u,B2u),B3u),conv(conv(A1u,A2u),conv(A3u,A4u)),500); %g*B1*B2*B3 / A1*A2*A3*A4
[H5u,W] = freqz(g_up*conv(conv(B1u,B2u),conv(B3u,B4u)),conv(conv(A1u,A2u),conv(A3u,A4u)),500); %g*B1*B2*B3*B4 / A1*A2*A3*A4


figure
plot(W,20*log10(abs(H1u)))
hold on
plot(W,20*log10(abs(H2u)))
hold on
plot(W,20*log10(abs(H3u)))
hold on
plot(W,20*log10(abs(H4u)))
hold on
plot(W,20*log10(abs(H5u)))

legend('H1','H2','H3','H4','H5')
title('Up Case: Magnitude Plots of Successive Cumulative Transfer Functions')
xlabel('Frequency')
ylabel('Magnitude (dB)')


% DOWN
B1d = sos_dwn(1,1:3);
B2d = sos_dwn(2,1:3);
B3d = sos_dwn(3,1:3);
B4d = sos_dwn(4,1:3);
A1d = sos_dwn(1,4:6);
A2d = sos_dwn(2,4:6);
A3d = sos_dwn(3,4:6);
A4d = sos_dwn(4,4:6);

[H1d,W] = freqz(g_dwn,A1d,500); %g / A1
[H2d,W] = freqz(g_dwn*B1d,conv(A1d,A2d),500); %g*B1 / A1*A2
[H3d,W] = freqz(g_dwn*conv(B1d,B2d),conv(conv(A1d,A2d),A3d),500); %g*B1*B2 / A1*A2*A3
[H4d,W] = freqz(g_dwn*conv(conv(B1d,B2d),B3d),conv(conv(A1d,A2d),conv(A3d,A4d)),500); %g*B1*B2*B3 / A1*A2*A3*A4
[H5d,W] = freqz(g_dwn*conv(conv(B1d,B2d),conv(B3d,B4d)),conv(conv(A1d,A2d),conv(A3d,A4d)),500); %g*B1*B2*B3*B4 / A1*A2*A3*A4

figure
plot(W,20*log10(abs(H1d)))
hold on
plot(W,20*log10(abs(H2d)))
hold on
plot(W,20*log10(abs(H3d)))
hold on
plot(W,20*log10(abs(H4d)))
hold on
plot(W,20*log10(abs(H5d)))

legend('H1','H2','H3','H4','H5')
title('Down Case: Magnitude Plots of Successive Cumulative Transfer Functions')
xlabel('Frequency')
ylabel('Magnitude (dB)')

%% (e) DO THIS PART!!!
%if (A1d == conv(conv(A1u,A2u),conv(A3u,A4u)))
 %   disp('Denominator polynomials appear in reverse order')
 
%% 5.
% Original Coefficient Values
H0 = [0.1336 0.0568 0.0563 0.1336; 1 -1.5055 1.2630 -0.3778];
HA = [-0.4954 1 0.7626 -1.0101 1; 1 -0.4954 1 -1.0101 0.7626];
b = 5; %number of total bits available

% y = fi(x,s,w,f) where x=value, s = 0 for unsigned and 1 for signed, w = wordlength, f = fraction length
% x0 = y.data returns numerical data
% Varying number of fractional bits
H01 = fi(H0,1,5,1); %really bad
H02 = fi(H0,1,5,2); %underflow for 2 values
H03 = fi(H0,1,5,3); %underflow for 2 values
H04 = fi(H0,1,5,4); %BEST
H05 = fi(H0,1,5,5); %identical to 4 but worse in 3 values

H0Q = H04;
x0 = H0Q.data
x0(2,1) = 1; %Pure 1's are not subject to the rounding as per the problem statement

HA1 = fi(HA,1,5,1); %all rounded to 1 or 0.5
HA2 = fi(HA,1,5,2); %all rounded to 1,0.75 or 0.5
HA3 = fi(HA,1,5,3); %same as HA2
HA4 = fi(HA,1,5,4); %same as HA2 but 1s are rounded (which we dont care about)
HA5 = fi(HA,1,5,5); %really bad

HAQ = HA3; %can select either HA2 or HA3 because they are identical
xA = HAQ.data

%%
%(B)
% H0
% Check ideal gain
H0_w0 = sum(H0(1,:))./sum(H0(2,:)) %analogous to plugging in z=1
H0_wpi = (H0(1,1) - H0(1,2) - H0(1,3) - H0(1,4))./(H0(2,1) - H0(2,2) - H0(2,3) - H0(2,4)) %analogous to plugging in z=-1 
% w=0 has an ideal gain of 1 and w=pi has an ideal gain of 0. There is
% slight rounding error present.

% Check quantized gain 
w = [0 pi]; %frequencies of interest
H0Q_freqdom = freqz(x0(1,:),x0(2,:),w); %num coeff,denom coeff
H0Q_magresp = abs(H0Q_freqdom) %Magnitude at w=0 and w=pi

% Report Error
H0Q_mag_dif = abs((1 - H0Q_magresp(1))); %absolute error - dont need to worry about w=pi portion because that is exact
H0Q_mag_error = 20*log10(H0Q_magresp(1)/1) %dB value of error

%%
% HA
% Check ideal gain
HA_w0 = 0.5*((HA(1,1) + 1)./(1 + HA(2,2)) + (HA(1,3) + HA(1,4) + 1)./(1 + HA(2,4) + HA(2,5))) %analogous to plugging in z=1
HA_wpi = 0.5*((HA(1,1) - 1)./(1 - HA(2,2)) + (HA(1,3) - HA(1,4) - 1)./(1 - HA(2,4) - HA(2,5))) %analogous to plugging in z=-1 
% w=0 has an ideal gain of 1 and w=pi has an ideal gain of 0. There is
% slight rounding error present.

% Check quantized gain 
w = [0 pi]; %frequencies of interest
HAQ1_freqdom = freqz(xA(1,1:2),2.*xA(2,1:2),w); 
HAQ2_freqdom = freqz(xA(1,3:5),2.*xA(2,3:5),w);
HAQ_freqdom = HAQ1_freqdom + HAQ2_freqdom;
H0Q_magresp = abs(HAQ_freqdom); %Magnitude at w=0 and w=pi

% Report Error
H0Q_mag_error = 20*log10(H0Q_magresp(1)/1) %dB value of error
% There is NO ERROR!! The values are exact :)

%% (c)
w_c = linspace(0, pi,1e4); %10^4 frequencies from 0 to pi

H_c = freqz(H0(1,:),H0(2,:),w_c); % Frequency Response for H(w) not quantized
Hq0_c = freqz(x0(1,:),x0(2,:),w_c); % Frequency Response for H0(w) quantized

HqA1_c = freqz(xA(1,1:2),2.*xA(2,1:2),w_c); 
HqA2_c = freqz(xA(1,3:5),2.*xA(2,3:5),w_c);
HqA_c = HqA1_c + HqA2_c; %Frequency Response for HA(w) quantized

% Maximum Values
MaxQ0 = max(abs(H_c - Hq0_c))
MaxQA = max(abs(H_c - HqA_c))

% Plot magnitude response of the three filters 
figure
plot(w_c, 20*log10(abs(H_c)))
hold on
plot(w_c, 20*log10(abs(Hq0_c)))
hold on
plot(w_c, 20*log10(abs(HqA_c)))
legend('H(w) not quantized', 'Hq0(w) original quantized', 'HqA(w) sum of two allpass quantized')
xlabel('frequency')
ylabel('Magnitude (dB)')
title('Magnitude Response of Three Filters')
ylim([-40 1])

%% (d)
%% 1 Maximum Deviation of the Filter Gains in Passband
% The infinite precision filter and the HqA filter are almost identical in
% the passband, preserving the amplitude response by keeping it around 0 dB
% with only a 0.92 dB deviation. The Hq0 filter has severe attenuation in
% the entire passband. At the best portion it attenuates by approximatly
% 3.5 dB and at the worst portion it attenuates by approximatly 10 dB. This
% is a total distortion of the passband amplitude and will ruin the signal.
passband_index = (1e4./pi).*(.3*pi);
H_pass = H_c(1:passband_index); % Passband of H(w) not quantized
Hq0_pass = Hq0_c(1:passband_index); % Passband of H0(w) quantized
HqA_pass = HqA_c(1:passband_index); % Passband of HA(w) quantized

Dev_Hq0 = 20*log10(abs(H_pass)) - 20*log10(abs(Hq0_pass));
Max_Dev_Hq0 = max(Dev_Hq0)

Dev_HqA = 20*log10(abs(H_pass)) - 20*log10(abs(HqA_pass));
Max_Dev_HqA = max(Dev_HqA)
% These computations show that quantitativly there is at most a 0.035
% difference between the infinite precision and parallel allpass filter
% while there is up to a 9.4947 (very large and significant) difference
% between the infinite precision filter and the directly quantized version.

%% 2 Deviation from the equiripple characteristics
% In passband, both the infinite precision filter and the HqA filter
% observe the equiripple requirement, while the Hq0 filter does not obey
% equiripple, with an amplitude monotonically attenuating. Therefore, there
% is a tolerance of approximately 6.5 dB in the stopband, which is
% qualitativly bad.

%% 3 Maximum Gain in the Stopband
% The maximum gain for the infinite precision filter is 20 dB. The HqA has
% a larger attenuation, which means it meets the specification. The Hq0 has
% smaller attenuation and does not meet the stopband specification.
stopband_index = 4139; %determined approximately from the graph
H_stop = H_c(stopband_index:end); % Stopband of H(w) not quantized
Hq0_stop = Hq0_c(stopband_index:end); % Stopband of H0(w) quantized
HqA_stop = HqA_c(stopband_index:end); % Stopband of HA(w) quantized

Max_Stopgain_H = max(20*log10(abs(H_stop)))
Max_Stopgain_Hq0 = max(20*log10(abs(Hq0_stop)))
Max_Stopgain_HqA = max(20*log10(abs(HqA_stop)))
% Notice the maximum gain in the stopband is approximately the same for the infinite
% precision with a value of -20/-21 dB. The maximum gain for directly
% quantized filter is -16dB which is much less attenuation.

%% All of these properties highlight the fact that using the sum of two allpass
% filters when quantizing will preserve the frequency response
% characteristics to a high degree of precision while quantizing the actual
% transfer function directly destroys the signal.

%% (e)
w_gd = linspace(0,0.3*pi,512); %passband only

H_gd = grpdelay(H0(1,:),H0(2,:),w_gd); %Group Delay for H(w) not quantized
Hq0_gd = grpdelay(x0(1,:),x0(2,:),w_gd); %Group Dealy for H0(w) quantized

HqA1_gd = grpdelay(xA(1,1:2),2.*xA(2,1:2),w_gd); 
HqA2_gd = grpdelay(xA(1,3:5),2.*xA(2,3:5),w_gd);
HqA_gd = HqA1_gd + HqA2_gd; %Group Dealy for HA(w) quantized

% Plot group delay of the three filters 
figure
plot(w_gd, H_gd)
hold on
plot(w_gd, Hq0_gd)
hold on
plot(w_gd, HqA_gd)
legend('H(w) not quantized', 'Hq0(w) original quantized', 'HqA(w) sum of two allpass quantized')
xlabel('frequency(rad)')
ylabel('Group Delay (dB)')
title('Group Delay of Three Filters')

% The parallel allpass filter performs worse when it comes to group delay.
% It has a larger group delay which is a very nonlinear function of
% frequency. It looks like an exponential curve. The less linear the group
% delay is (which does NOT mean this is a non linear FILTER), the more 
% phase distortion in the time domain.

%% (f)
[z_H,p_H,k_H] = tf2zpk(H0(1,:),H0(2,:)) %Zeros, Poles for H(w) not quantized
[z_H0q,p_H0q,k_H0q] = tf2zpk(x0(1,:),x0(2,:)) %Zeros, Poles for H0(w) quantized

tf1 = tf(xA(1,1:2),2.*xA(2,1:2),0.1,'Variable','z^-1');
tf2 = tf(xA(1,3:5),2.*xA(2,3:5),0.1,'Variable','z^-1');
tf_HqA = tf1 + tf2;
[z_HAq,p_HAq,k_HAq] = tf2zpk(cell2mat(tf_HqA.num),cell2mat(tf_HqA.den));

figure
subplot(1,3,1)
zplane(z_H,p_H)
title('Poles and Zeros of Unquantized H(w)')

subplot(1,3,2)
zplane(z_H0q,p_H0q)
title('Poles and Zeros of Quantized H0(w)')

subplot(1,3,3)
zplane(z_HAq,p_HAq)
title('Poles and Zeros of Quantized HA(w) - Parallel AllPass')