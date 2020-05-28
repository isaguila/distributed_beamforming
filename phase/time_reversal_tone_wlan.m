close all;
clear;
clc;

rng(0);

set(0, 'DefaultFigureWindowStyle', 'docked');

snr = 25;

wlan_chan_model = wlanTGnChannel(...
    'PathGainsOutputPort', true, ...
    'DelayProfile', 'Model-B');

fs = wlan_chan_model.SampleRate;
fc = fs/16;
tone_len = 2^9;
t = 0:tone_len-1;
t = t(:)/fs;
x = exp(1i*2*pi*fc*t);

t = 1:tone_len;
t = t(:);

[~, path_gains] = wlan_chan_model(x);

h = path_gains(1, :).';
phase_response = angle(h);
mag_response = 20*log10(abs(h));

y = conv(x, h, 'same');
y = y*norm(x)/norm(y);

x_tr = conv(x, conj((flip(h))), 'same');
y_tr = conv(x_tr, h, 'same');
y_tr = y_tr*norm(x)/norm(y_tr);

p_noise = mean(abs(y_tr).^2)/(10^(snr/10));
noise1 = sqrt(p_noise/2)*randn(length(y), 2)*[1; 1i];
noise2 = sqrt(p_noise/2)*randn(length(y_tr), 2)*[1; 1i];
y = y + noise1;
y_tr = y_tr + noise2;

y = y*norm(x)/norm(y);
y_tr = y_tr*norm(x)/norm(y_tr);

figure(1);
ax1 = subplot(3, 1, 1);
plot(t, real(x), t, imag(x), 'LineWidth', 2);
legend('real', 'imag');
xlabel('Sample Index (n)');
ylabel('Amplitude');
title('x(t)');
grid on;

figure(1);
ax2 = subplot(3, 1, 2);
plot(t, real(y), t, imag(y), 'LineWidth', 2);
legend('real', 'imag');
xlabel('Sample Index (n)');
ylabel('Amplitude');
title('y(t)');
grid on;

figure(1);
ax3 = subplot(3, 1, 3);
plot(t, real(y_tr), t, imag(y_tr), 'LineWidth', 2);
legend('real', 'imag');
xlabel('Sample Index (n)');
ylabel('Amplitude');
title('y_{tr}(t)');
grid on;

linkaxes([ax1, ax2, ax3], 'xy');