%% Shifra Abittan
% DSP HW #6: FIR Filter Design

%% 1(a)
ChebWin = chebwin(31,30);

figure
t = 0:1:30;
stem(t,ChebWin)
title('31 Point CHEBYCHEV WINDOW with 30dB peak sidelobe level in TIME DOMAIN')
xlabel('time')
%Already normalized?

%% 1(b)
zplane(ChebWin.')

%% 1(c)
ChebWinFreq = fft(ChebWin,1000);
ChebWinFreq = fftshift(ChebWinFreq);
Fwin = linspace(-1,1,1000);
figure 
plot(Fwin,20*log10(abs(ChebWinFreq)))
title('Magniutde Response of a 31 Point Chebychev Window')
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')

%% 1(d)

% Generate Kaiser Window
beta = 3.14; %3.15,3.16 (both values line up nulls perfectly) look at zero values
KaiserWin = kaiser(31,beta);

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
KaiserFreq = fft(KaiserWin./31.6228,1000);
KaiserFreq2 = fftshift(KaiserFreq);

figure 
plot(Fwin,20*log10(abs(ChebWinFreq)))
hold on
plot(Fwin,20*log10(abs(KaiserFreq2)))
title('Frequency Response of Window Functions with Same Main Lobe Widths')
legend('Chebychev Window','Kaiser Window')

figure
zplane(ChebWin.',KaiserWin.')

%% 1(g)
%PeakSideLobeLevel = min(20*log10(W(0)/mag(W(w))))
% min from first null to pi

%% 1(h) Fraction of Energy in Sidelobes
% Total Energy = mag squared added up of w[n]
% Fraction = sum of mag sq of freq coefficients from first null to pi
% Chebychev Window
TotalE_Cheb = sum((abs(ChebWin)).^2);
%SidelobeE_Cheb = sum((abs(ChebWinFreq(:))).^2) %insert index from null
%FractionCheb = SidelobeE_Cheb./TotalE_Cheb

TotalE_Kaiser = sum((abs(KaiserWin)).^2);
%SidelobeE_Kaiser = sum((abs(KaiserFreq(:))).^2) %insert index from null
%FractionKaiser = SidelobeE_Kaiser./TotalE_Kaiser

%% 2(b)
PDev = (10.^(2/20) - 1)./(10^(2/20) + 1);
SDev = 10^(-30/20);
%formulas may be wrong!

%% 2(c)
% KAISER DESIGN
F = [1.5*10e6 2*10e6 3*10e6 3.5*10e6]; %band edge freq in Hz
A = [0 1 0]; %ideal behavior in band
DEV = [SDev,PDev,SDev]; %maximum deviations on linear scale
Fs = 10*10e6; %sampling frequency

[N,Wn,BTA,FILTYPE] = kaiserord(F,A,DEV,Fs);
B = fir1(N,Wn,FILTYPE,kaiser(N+1,BTA),'noscale');

% PARKS-MCCLELLAN EQUIRIPPLE DESIGN
[N_rip,Fo_rip,Ao_rip,W_rip] = firpmord(F,A,DEV,Fs);
B_rip = firpm(N_rip,Fo_rip,Ao_rip,W_rip);

%% 2(d)
% Specify Filter Orders
KaiserOrder = N
EquirippleOrder = N_rip

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

%% 2(e)
InvProp = W_rip.*DEV.'
% This statement above shows that the weights are inversely proportional to
% the deviations in the bands. Mathematically W prop to 1/deviations.
% Therefore, W.*deviations should equal a constant. InvProp is equal to a
% constant value for all the weight/deviation pairs.