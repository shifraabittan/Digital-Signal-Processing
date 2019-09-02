% Shifra Abittan
% DSP HW4: DFT

%% (2)
% (A)
N = 256; %DFT block size
fs = 50e3; %sampling freq
binspacing = fs/N %in Hz

% (B)
k1 = 10e3*N/fs; 
k2 = 40e3*N/fs;

k = [round(k1) round(k2)]

% (C) Straddle Loss

inputfreq = [0.2*2*pi/N -0.2*2*pi/N]; %use radian binfrequency
% Use 0.2/-0.2 because this is the distance between actual location of
% 10kHz/40kHz delta and the closest bin frequency

value = diric(inputfreq,256); %Windowing is generally an independent stage from
% the DFT but the diriclet delta here represents the DFT output because the
% output is just a sin wave multiplied with a rectangular truncation window = 
% convolution in frequency domain of delta at sin freq with a diriclet sinc =
% a diriclet sinc shifted to the location of the delta/frequency of the sin
% wave.
straddleloss = 20*log10(value/1) %in dB

% (D) Hamming Window
window = hamming(256);
W_0 = abs(sum(window)) %using DTFT formula directly. The exponential term 
% drops out because e^0 = 1 and w = 0 here

% compute W(w = 0.2*binspacing) using DTFT formula
w = 0.2*2*pi/256;

index = 0:1:255;
exponential = exp(-j.*w.*index); %compute exponential term corresponding to each x[n]
timeterms = window.*exponential.';

W_w = abs(sum(timeterms))

ham_straddleloss = 20*log10(W_w/W_0)

% (E) N=512
E_N = 512;
E_inputfreq = [0.2*2*pi/E_N -0.2*2*pi/E_N]; %use radian binfrequency
% Use 0.2/-0.2 because this is the distance between actual location of
% 10kHz/40kHz delta and the closest bin frequency

E_value = diric(E_inputfreq,512);
E_straddleloss = 20*log10(E_value/1) %in dB

% Hamming Window
E_window = hamming(512);
E_W_0 = abs(sum(E_window)); %using DTFT formula directly. The exponential term 
% drops out because e^0 = 1 and w = 0 here

% compute W(w = 0.2*binspacing) using DTFT formula
E_w = 0.2*2*pi/512;

E_index = 0:1:511;
E_exponential = exp(-j.*E_w.*E_index); %compute exponential term corresponding to each x[n]
E_timeterms = E_window.*E_exponential.';

E_W_w = abs(sum(E_timeterms));

E_ham_straddleloss = 20*log10(E_W_w/E_W_0)

%% (6)
fs = 100e6;
T = 1/fs;

t = 0:T:999*T;
sinwave = 2*sin(2*pi*20e6*t);
awgn = sqrt(0.2)*randn(1,1000);
corrupted = sinwave + awgn;

Cheb_Window = chebwin(1000,30);

windowed = corrupted.*Cheb_Window.';

DFT = fft(windowed,1024);
DFT2= fftshift(DFT);
k = -1024/2:1024/2-1;
w = 2.*pi.*k/1024;
f = fs.*k/1024;

figure
plot(f,20*log10(abs(DFT2)))
title('DFT Coefficient Magnitudes for a Noisy Sinwave with a Chebychev window')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%% (7)
% (A)
w_0 = repmat([0 0 0 1],1, 250);
applywindow = corrupted.*w_0;

DFT7 = fft(applywindow,1024);
DFT7= fftshift(DFT7);
k = -1024/2:1024/2-1;
w = 2.*pi.*k/1024;
f = fs.*k/1024;

figure
plot(f,20*log10(abs(DFT7)))
title('Spectrum of Noisy Sinwave with a w_0 window')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

% There are more peaks (approximatly 4 for every 1 peak before). The peaks
% appear in a periodic form with a distance of approximatly 1 then 1.5 then 
% 1 etc between them. The frequencies where they appear are: 
% -4.5, -3, -2, -0.5, 2, 3, 4.5 Hz.


% (B)
DFT_w0 = fft(w_0,1000);
k = 0:999;
figure 
stem(k, abs(DFT_w0))
title('Stem Plot of the Magnitude of 1000 Point DFT on w_0[n]')
xlabel('k (DFT index)')
ylabel('Magnitude (linear)')
% Four deltas can be seen. This means all frequencies are
% attenuated, with no energy, except these four points.
% In the spectrum, there were four tall delta peaks, which correspond to the
% four peaks seen here. This makes sense because the Fourier Transform of a 
% deltra train is a scaled and shifted delta train. This is why the deltas 
% have different widths between them and different heights.
