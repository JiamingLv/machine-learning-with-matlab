%%
clear
clc
close all
% % %
% % Fs = 1000;                    % Sampling frequency
% % T = 1/Fs;                     % Sampling period
% % L = 1024;                     % Length of signal
% % t = (0:L-1)*T;                % Time vector
% % x1 = cos(2*pi*50*t);          % First row wave
% % plot(t, x1)
% % r = 2^nextpow2(L);
% % % DFT
% % -----------自编DFT-------------%
% % sum = 0;
% % N = L;
% % for k = 0:N-1
% %     for n = 0:N-1
% %         W = exp(-1i * 2*pi / N);
% %         sum = sum + x1(n+1) * W^(k*n);
% %     end
% %     y(k+1) = sum;
% %     sum = 0;
% % end
% % -------------------------------%
% % y1 = fft(x1, r);
% % % 显示
% % f = Fs * (0 : (r-1)) / r;
% % P = abs(y/N*2);
% % P1 = abs(y1/N*2);
% % figure
% % subplot(2,1,1)
% % plot(f, P(1:r))
% % subplot(2,1,2)
% % plot(f, P1(1:r/2+1))

%%

Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sampling period
L = 1024;                     % Length of signal
t = (0:L-1)*T;                % Time vector
x1 = cos(2*pi*50*t);          % First row wave

Nt = L;
times = log2(Nt);
index = 0:Nt-1;
index_y = zeros(Nt,1);
%------------------自编FFT-------------------%
%倒序
for j = 1:Nt
    for i = 1:times
        index_y(j) = bitshift(index_y(j),1); %y左移一位
        index_y(j) = index_y(j) + rem(index(j),2); %y加上x对2的余数
        index(j) = bitshift(index(j),-1); %x右移一位
    end
end

%FFT
x = zeros(Nt,1);
x(1:Nt) = x1(index_y(1:end)+1);
xs = x;
for i = 1:times
    m = 2^(i-1);
    for j = 1:2^(i):Nt
        for k = j:j+m-1
            ts = xs(k+m) * exp(-1i * 2*pi * (k-j)/2^(i));
            y(m+k) = xs(k) - ts;
            y(k) = xs(k) + ts;
        end
    end
    xs = y;
end
%------------------------------%
r = L;
y1 = fft(x1, r);
%% 显示
f = Fs * (0 : (r-1)) / r;
P = abs(y/L*2);
P1 = abs(y1/L*2);
figure
subplot(2,1,1)
plot(f, P(1:r))
subplot(2,1,2)
plot(f, P1(1:r))
