%% Shifra Abittan
% DSP HW #5: MultiDimensional Signals and Systems

%% 2
% (b)
h = (1/6)*[1 4 1; 4 -20 4; 1 4 1]; %2D FIR filter impulse response in time
[H,Fx,Fy] = freqz2(h,100); %frequency domain(dimensions to use)??

H_MaxImag = max(abs(imag(H)));

%Max Numerical Difference between freq2 and formula using 100x100 discrete
%points
w1 = linspace(-pi,pi,100);
w2 = linspace(-pi,pi,100);
[W1,W2] = meshgrid(w1,w2);

Hformula = -20/6 + (8/6)*cos(W1) + (8/6)*cos(W1) + (1/3)*cos(W1+W2) + (1/3)*cos(W1-W2);

NumError = max((max(abs(Hformula - real(H)))).') 

%Alternatively, because the w1,w2 does not hit every possible
%value(including pi,pi - which happens to be the maximum):
%Hformula(pi,pi) = 5.3 (absolute value)
MaxofH = max(max(abs((real(H)).')));
NumErrorAlt = 5+(1/3) - MaxofH;

H = real(H);

% (c)
% Contour Plot
figure 
sgtitle('Contour Plot of the Frequency Response H(w1,w2)')
contour(w1,w2,H)
xlabel('w1 (rad)')
ylabel('w2 (rad)')
% Axes Radian Markings in Terms of Pi
set(gca,'XTick',-pi:pi/4:pi) 
set(gca,'XTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',-pi:pi/4:pi) 
set(gca,'YTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
 
% Surface Plot
figure
sgtitle('Surface Plot of the Frequency Response H(w1,w2)')
surf(w1,w2,H)
xlabel('w1 (rad)')
ylabel('w2 (rad)')
% Axes Radian Markings in Terms of Pi
set(gca,'XTick',-pi:pi/4:pi) 
set(gca,'XTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',-pi:pi/4:pi) 
set(gca,'YTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})

%The filter looks like a lowpass filter. It allows frequencies close to 0
%through and attenuates at higher frequencies. Near zero, the frequency
%response is 0, meaning the original signal is untouched. 

%(d)
figure
laplacian = -(Fx.^2);
H_slice = H(50,:);

plot(Fx,H_slice);
hold on
plot(Fx,laplacian);
legend('H(wx,0):Slice','Ideal Curve/Laplacian')
xlabel('Frequency from -\pi to \pi, normalized from -1 to 1')

%(e)
[ff1,ff2] = meshgrid(Fx,Fx);
spec = -(ff1.^2 + ff2.^2)*pi;
figure
surf(Fx,Fy,spec)
hold on
surf(Fx,Fy,H)
colormap hsv

% Axes Radian Markings in Terms of Pi
set(gca,'XTick',-pi:pi/4:pi) 
set(gca,'XTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',-pi:pi/4:pi) 
set(gca,'YTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})

legend('Spectrum of Laplacian','Plot of H')


%% 3
%(a) see function upsamp below
%(b)
g = upsamp(h);
[Hup,Fxup,Fyup] = freqz2(g,100);
figure
subplot(1,2,1)
contour(Hup)
subplot(1,2,2)
contour3(Hup)

figure
subplot(1,2,1)
surface(Hup)
subplot(1,2,2)
surf(Hup)

% This is a plausible approximation; however, because of the sampling, there
% is periodization. Lowpass filtering would select just the filter in the
% middle, which resembles that of part c.


%% 4
%(a)
%Impulse Time Domain
hy = (1/8)*[-1 -2 -1; 0 0 0; 1 2 1];
hx = hy.';

%Frequency Domain Hy
[Hy,Fx,Fy] = freqz2(hy,100);
figure
surf(abs(Hy))
title('Surface Plot of the Frequency Response Hy')

%Frequency Domain Hx
[Hx,Fx,Fy] = freqz2(hx,100);
figure 
surf(abs(Hx))
title('Surface Plot of the Frequency Response Hy')

% These are vertical and horizontal highpass

%(c)
cirimag = imread('circuit.tif');
convdouble = im2double(cirimag);

image(cirimag)
EdgesL1 = SobolEdge(convdouble,0.0555,hy,hx,1); %1=L1,2=L2
imtool(EdgesL1)

medvalue = median(convdouble(:))
%%
%(d)
EdgesL2 = SobolEdge(convdouble,0.0555,hy,hx,2);
imtool(EdgesL2)

%%
function y = upsamp(x)
    [r,c] = size(x);
    
    % Insert Zeros into Columns
    temp = [x;zeros(size(x))];
    emptycol = reshape(temp,r,c*2);
    emptycol = emptycol(:,1:end-1); %remove the last row of zeros
    
   
    % Insert Zeros into Rows
    [r1,c1] = size(emptycol);
    temp = [emptycol zeros(size(emptycol))];
    tempT = temp.';
    emptyrow = reshape(tempT,c1,r1*2);
    emptyrowT = emptyrow.';
    y = emptyrowT(1:end-1,:);

end

%% 4.b Compute the Edges based on the Approximatation of Gradient
function Edges = SobolEdge(g,T,hy,hx,L)
    %Compute Convolutions of signal with hy and hx to approximate the
    %partial derivatives
    del_y = conv2(hy,g);
    del_x = conv2(hx,g);
    
    % There is 1 extra row/column on each side. This is because the size
    % of the result of convolution is N+M-1 where N is the size of matrix 1
    % and M is the size of matrix 2. Because hy and hx are 3x3, N+3-1 = N+2 
    % in the row and column direction. This is equivalent to 1 row/column
    % added on each of the four sides.
    [ry,cy] = size(del_y);
    del_y = del_y(2:ry-1,2:cy-1);
    
    [rx,cx] = size(del_x);
    del_x = del_x(2:rx-1,2:cx-1);
    
    if L == 1
        contrast = ((del_y.^2 + del_x.^2).^(1/2));
    end
    
    if L == 2
        contrast = (del_y + del_x)
    end
    Edges = contrast>T;
end
    
    