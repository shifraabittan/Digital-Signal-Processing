 %% Shifra Abittan
% Professor Fontaine
% DSP ECE310
% Problem Set 1: Transfer Function Analysis

% (2)
%% Check to ensure proper all-pass and minimum-phase transfer functions were found
figure 
subplot(2,2,1)
H_z = [-4/3; -5];
H_p = [1/2; -3; -3];
zplane(H_z, H_p)
title('Zero-Pole Plot of H(z)')

subplot(2,2,2)
A_z = [-3/4; -1/5; -3; -3; 0];
A_p = [-4/3; -5; -1/3; -1/3];
zplane(A_z, A_p)
title('Zero-Pole Plot of A(z)')

subplot(2,2,3)
HA_z = [H_z; A_z];
HA_p = [H_p; A_p];
zplane(HA_z, HA_p)
title('Zero-Pole Plot of H(z)*A(z)')

subplot(2,2,4)
Hmin_z = [-3/4; -1/5; 0];
Hmin_p = [1/2; -1/3; -1/3];
zplane(Hmin_z, Hmin_p)
title('Zero-Pole Plot of Hmin(z)')

%% Graph three phase responses

% H(z)
[b_H,a_H] = zp2tf(H_z,H_p,1);
[H,W] = freqz(b_H,a_H,1000);

% A(z)
[b_A,a_A] = zp2tf(A_z,A_p,1);
[H_A,W_A] = freqz(b_A,a_A,1000);

% Hmin(z)
[b_Hmin,a_Hmin] = zp2tf(Hmin_z,Hmin_p,1);
[H_Hmin,W_Hmin] = freqz(b_Hmin,a_Hmin,1000);

figure
plot(W,unwrap(angle(H))*180/pi,W_A,unwrap(angle(H_A))*180/pi,W_Hmin,unwrap(angle(H_Hmin))*180/pi)
legend('H(z)', 'A(z)', 'Hmin(z)')
title('Phase Responses')
xlabel('Normalized Frequency') 
ylabel('Unwrapped Phase Response (degrees)')
%Note that Hmin has "less phase" than H as expected

%% (4)
% Part A Result: S(z) = (-2*z^2 + 5z - 2)^2 / z*(6*z^2 + 13*z + 6)
% Spectral Factorization: (1) Find poles and zeros
syms z;

% zeros
S_num = (-2*z^2 + 5*z - 2)^2;
S_num = expand(S_num);
zeros = roots(coeffs(S_num))

%poles
S_denom = z*(6*z^2 + 13*z + 6);
S_denom = expand(S_denom);
poles = roots(coeffs(S_denom)) %dont forget the pole at 0

%figure
%zplane(coeffs(S_num), coeffs(S_denom))

% Spectral Factorization: (2) Innovations filter H(z) = poles/zeros in
% |z|<1. Whitening filter G(z) = remaining poles/zeros = paraconjugate of
% H(z)

zeros_H = [1/2 1/2]
poles_H = [-2/3 0]





