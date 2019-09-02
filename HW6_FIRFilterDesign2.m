%% Shifra Abittan
% DSP HW #6: FIR Filter Design

%% 1(a)
ChebWin = chebwin(31,30);
ChebWin = ChebWin./sum(ChebWin); %normalize

figure
t = 0:1:30;
stem(t,ChebWin)
title('31 Point CHEBYCHEV WINDOW with 30dB peak sidelobe level in TIME DOMAIN')
xlabel('time')

%% 1(b)
zplane(ChebWin.')

[ChebWinFreq,W] = freqz(ChebWin,1,1000,2);

%Identify first null
Null = (min(find(ChebWinFreq<0.02)))
MainLobeWidth = 2*(Null/1000)/(4*pi/31)

%% 1(c)
%ChebWinFreq2 = fft(ChebWin/sum(ChebWin),1000);
%ChebWinFreq = fftshift(ChebWinFreq);
%Fwin = linspace(-1,1,1000);

figure 
f = linspace(0,1,1000); %frequency on a normalized scale from 0 to 1/Nyquist bandwidth
plot(f,20*log10(abs(ChebWinFreq)))
title('Magnitude Response of a 31 Point Chebychev Window')
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')

%% 1(d)

% Generate Kaiser Window
beta = 3.14; 
KaiserWin = kaiser(31,beta);
KaiserWin = KaiserWin./sum(KaiserWin);

% Plot Time Domain
%{
figure
t = 0:1:30;
stem(t,KaiserWin)
title('31 Point KAISER WINDOW with b= in TIME DOMAIN')
xlabel('time')
%}

%% 1(e)
% Superimpose Freq Response of Kaiser on Chebychev Window to find beta

%KaiserFreq = fft(KaiserWin./31.6228,1000);
%KaiserFreq2 = fftshift(KaiserFreq);

KaiserFreq = freqz(KaiserWin,1,1000,2);

figure 
plot(f,20*log10(abs(ChebWinFreq)))
hold on
plot(f,20*log10(abs(KaiserFreq)))
title('Frequency Response of Window Functions with Same Main Lobe Widths')
legend('Chebychev Window','Kaiser Window')

figure
zplane(ChebWin.',KaiserWin.')

%% 1(f)
beta

%% 1(g)
PeakSideLobeLevel = min(20*log10(KaiserFreq(1)/min(abs(KaiserFreq(34:end)))))


%% 1(h) Fraction of Energy in Sidelobes
% Total Energy = mag squared added up of w[n]
% Fraction = sum of mag sq of freq coefficients from first null to pi
% Chebychev Window
TotalE_Cheb = sum((abs(ChebWin)).^2);
SidelobeE_Cheb = sum((abs(ChebWinFreq(34:end))).^2); %insert index from null
FractionCheb = SidelobeE_Cheb./TotalE_Cheb

Kaisernull = min(find(KaiserFreq<0.02)); %this is index 34
TotalE_Kaiser = sum((abs(KaiserWin)).^2);
SidelobeE_Kaiser = sum((abs(KaiserFreq(34:end))).^2); %insert index from null
FractionKaiser = SidelobeE_Kaiser./TotalE_Kaiser

%% 2(b)
PDev = (10.^(2/20) - 1)./(10^(2/20) + 1);
SDev = 10^(-30/20);

%% 2(c)
% KAISER DESIGN
F = [1.5e6 2e6 3e6 3.5e6]; %band edge freq in Hz
A = [0 1 0]; %ideal behavior in band
DEV = [SDev,PDev,SDev]; %maximum deviations on linear scale
Fs = 10e6; %sampling frequency

[N,Wn,BTA,FILTYPE] = kaiserord(F,A,DEV,Fs);
B = fir1(N+1,Wn,FILTYPE,kaiser(N+2,BTA),'noscale');

% PARKS-MCCLELLAN EQUIRIPPLE DESIGN
[N_rip,Fo_rip,Ao_rip,W_rip] = firpmord(F,A,DEV,Fs);
B_rip = firpm(N_rip+4,Fo_rip,Ao_rip,W_rip); 

%% 2(d)
% Specify Filter Orders
KaiserOrder = N+1
EquirippleOrder = N_rip+4

% Time Stem Plots

% Pad shorter filter with zeros (generic code that allows for any
% difference in length between the two filters) 
[~,c] = size(B);
[~,c_rip] = size(B_rip);

pad = abs(c-c_rip);
% check which filter needs padding
if c-c_rip < 0
    B = [B zeros(1,pad)];
end
if c-c_rip > 0
    B_rip= [B_rip zeros(1,pad)];
end

tfilt = 0:c-1;
figure
subplot(1,2,1)
stem(tfilt,B)
title('Kaiser Digital Bandpass Filter in Time Domain')
xlabel('time')
xlim([0 c-1])

subplot(1,2,2)
stem(tfilt,B_rip)
title('Parks-McClellan Equiripple Digital Bandpass Filter in Time Domain')
xlabel('time')
xlim([0 c-1])

% Magnitude Response in Frequency

[H,~] = freqz(B);
[H_rip,~] = freqz(B_rip);

F = linspace(0,5*10e6,512);

figure
subplot(1,2,1)
plot(F,20*log10(abs(H)))
title('Magnitude Response of Kaiser Filter')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-50 2])

% Plot Specs
hold on
plot(repmat(2*10e6,1,53),(-50:1:2),'--','Color','r') %passband
hold on
plot(repmat(3*10e6,1,53),(-50:1:2),'--','Color','r') %passband
hold on
plot(repmat(1.5*10e6,1,53),(-50:1:2),'--','Color','g') %stopband
hold on
plot(repmat(3.5*10e6,1,53),(-50:1:2),'--','Color','g') %stopband
hold on
plot(F,repmat(0,1,512),':','Color','k') %passband tolerance
hold on
plot(F,repmat(-2,1,512),':','Color','k') %passband tolerance

subplot(1,2,2)
plot(F,20*log10(abs(H_rip)))
title('Magnitude Response of Parks-McClellan Equiripple Filter')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-50 2])

% Plot Specs
hold on
plot(repmat(2*10e6,1,53),(-50:1:2),'--','Color','r') %passband
hold on
plot(repmat(3*10e6,1,53),(-50:1:2),'--','Color','r') %passband
hold on
plot(repmat(1.5*10e6,1,53),(-50:1:2),'--','Color','g') %stopband
hold on
plot(repmat(3.5*10e6,1,53),(-50:1:2),'--','Color','g') %stopband
hold on
plot(F,repmat(0,1,512),':','Color','k') %passband tolerance
hold on
plot(F,repmat(-2,1,512),':','Color','k') %passband tolerance

%In both the Kaiser design and the PM Equiripple design, the stopband 
%attentuation matches the specifications of 30dB down. The passband
%specification is not met perfectly in either but it is closer than in
%lower orders. For the PM Equiripple, any order lower than 23 does not meet
%the stopband specs. I tried increasing the order tremendously and it does
%not help with the passband specification failure in any specific way.
%% 2(e)
InvProp = W_rip.*DEV.'
% This statement above shows that the weights are inversely proportional to
% the deviations in the bands. Mathematically W prop to 1/deviations.
% Therefore, W.*deviations should equal a constant. InvProp is equal to a
% constant value for all the weight/deviation pairs.