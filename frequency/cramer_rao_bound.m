close all;
clear;
clc;

fsize = 14;

sigma = 1;  % noise variance (var(real) = var(imag) = sigma^2)
n0 = 1;
fs = 20e6;  % Based on WLAN
T = 1/fs;  % sampling period
b0 = 1; % amplitude
% snr = b0/2n0^2

snr_db = 0:50;
snr = 10.^(snr_db/10);

N = 500;
crb = 6./(snr*T^2*N.*(N^2-1));
semilogy(snr_db, crb, 'LineWidth', 2, 'DisplayName', ['N = ', num2str(N)]);
hold on;
xlabel('SNR (dB)', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Root Mean Square Error (Hz)', 'FontSize', fsize, 'Interpreter', 'Latex');
title('Cramer Rao Bound for One Shot Frequency Estimators', 'FontSize', fsize, 'Interpreter', 'Latex');
grid on;

N = 1000;
crb = 6./(snr*T^2*N.*(N^2-1));
semilogy(snr_db, crb, 'LineWidth', 2, 'DisplayName', ['N = ', num2str(N)]);
legend('show', 'FontSize', fsize);

N = 5000;
crb = 6./(snr*T^2*N.*(N^2-1));
semilogy(snr_db, crb, 'LineWidth', 2, 'DisplayName', ['N = ', num2str(N)]);
legend('show', 'FontSize', fsize);

N = 10000;
crb = 6./(snr*T^2*N.*(N^2-1));
semilogy(snr_db, crb, 'LineWidth', 2, 'DisplayName', ['N = ', num2str(N)]);
legend('show', 'FontSize', fsize);
set(gcf, 'color', 'w');

fn = '/Users/Ivan/Documents/Thesis/figures/crb_freq';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');

