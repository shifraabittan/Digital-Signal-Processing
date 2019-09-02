%% Shifra Abittan
% DSP HW #7: Multirate Systems

%% 4.
h = [1 2 3 4 5]; %Vector used to test
syms z
Hz = expand(poly2sym(h,z)/z^(length(h)-1));
% This line of code uses a vector h as the coefficients of a polynomial in
% z. The last h value corresponds to z^0 and the first value corresponds to
% z^(L-1) where L is the length of vector h. The polynomial is actually
% displayed with negative powers of z instead.

% Autonomous Function to Convert h into H(z)
Hz_fun = @(h) expand(poly2sym(h,z)/z^(length(h)-1));
Hz_fun(h) %check that autonomous function works

% Paraconjugate - for the z domain, the paraconjugate can be found by
% simply reversing the coefficients
Paraconj = @(Hz) expand(poly2sym(sym2poly(subs(Hz,z,z^-1)),z)/z^(length(sym2poly(subs(Hz,z,z^-1)))-1));
Paraconj(Hz)


%% 5.
% Coefficients of an FIR filter H0(z)
h = zeros(1,12);
h(1) = 0.15774243;
h(2) = 0.69950381;
h(3) = 1.06226376;
h(4) = 0.44583132;
h(5) = -0.31998660;
h(6) = -0.18351806;
h(7) = 0.13788809;
h(8) = 0.03892321;
h(9) = -0.04466375;
h(10) = -7.83251152e-4;
h(11) = 6.75606236e-3;
h(12) = -1.52353381e-3;

N=12; %order of FIR filters

H0 = h;

H1 = zeros(1,12);
F0 = zeros(1,12);
F1 = zeros(1,12);

for k = 0:11
    H1(k+1) = ((-1)^k)*h(N-k);
    F0(k+1) = h(N-k);
    F1(k+1) = ((-1)^k)*h(k+1);
end

%% (a)
[H0_f,W0] = freqz(H0,1);
[H1_f,W1] = freqz(H1,1);

figure
plot(W0,abs(H0_f),W1,abs(H1_f))
xlim([0 pi])
legend('|H0(w)|','|H1(w)|')
title('Magnitude Response of Two Filters H0 and H1')
xlabel('w (rad)')

%% (b)
[F0_f,~] = freqz(F0,1);
[F1_f,~] = freqz(F1,1);

figure
plot(W0,abs(F0_f),W1,abs(F1_f))
xlim([0 pi])
legend('|F0(w)|','|F1(w)|')
title('Magnitude Response of Two Filters F0 and F1')
xlabel('w (rad)')

% Confirm that mag of F(w) = H(w)
if sum(abs(H0_f)-abs(F0_f)) < 0.000001 %this threshhold is necessary because of roundoff in the original h coefficients
    disp('|F1(w)| = |H1(w)|')
end

if sum(abs(H1_f)-abs(F1_f)) < 0.000001
    disp('|F2(w)| = |H2(w)|')
end

%% (c)
Constant = abs(H0_f).^2 + abs(H1_f).^2;
% This is approximatly a constant ~ 4. Again, because of the roundoff in
% the original h coefficients, especially because I manually entered them
% instead of downloading from a library with higher precision, the constant
% is approximate.

%% (d)
if sum(abs(H1_f)-flip(abs(H0_f))) < 0.000001
    disp('H0 and H1 are quadrature mirror filters')
end

% This code takes H0 and effectivly shifts it by pi/2 and then ensures that
% the filters match. Because filters are periodic, a shift by pi/2 is the
% same as flipping the coefficients from top to bottom because an image
% shifts into the region from 0 to pi/2.

%% (e) 
%Time Domain Coefficients
e00 = H0(1:2:end); %even
e01 = H0(2:2:end); %odd
e10 = H1(1:2:end); %even
e11 = H1(2:2:end); %odd

%Z transform
E00 = Hz_fun(e00);
E01 = Hz_fun(e01);
E10 = Hz_fun(e10);
E11 = Hz_fun(e11);

E = [E00 E01;E10 E11];

%Paraconjugate
E00P = Paraconj(E00);
E01P = Paraconj(E01);
E10P = Paraconj(E10);
E11P = Paraconj(E11);

E_para = [E00P E10P;E01P E11P];%Hermitian Paraconjugate - each element 
% paraconjugated and the elements transposed

Paraunitary = E_para*E;
vpa(expand(simplify(Paraunitary)),2)

% If you examine the Paraunitary simplified results, the matrix is a diagonal
% matrix with values on the diagonal and zeros on the off diagonal. In the 
% non-zero valued entries, all of the terms are really tiny, effectivly zero, 
% except for one term in each which has the value 2*z^-5. Therefore, c=2
% and the delay element is 5.

%% (f)
%Time Domain Coefficients
r00 = F0(2:2:end); %odd F0
r01 = F1(2:2:end); %odd F1
r10 = F0(1:2:end); %even F0
r11 = F1(1:2:end); %even F1

%Z transform
R00 = Hz_fun(r00);
R01 = Hz_fun(r01);
R10 = Hz_fun(r10);
R11 = Hz_fun(r11);

R = [E00 E01;E10 E11];

PARTF = R*E;
vpa(expand(simplify(PARTF)),2)
%vpa(expand(simplify(vpa(expand(simplify(PARTF)),2)./2)),2)

