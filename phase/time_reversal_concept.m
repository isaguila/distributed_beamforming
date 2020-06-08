close all;
clear;
clc;

bits = randi([0 1], 10000, 1);
qam = qammod(bits, 4, 'UnitAveragePower', true, 'InputType', 'bit');

num_taps = 5;
h = randn(num_taps, 1);
h(1) = 1;
h(2:4) = 0;
h_tr = flip(h);

h_eq = conv(h, h_tr, 'same');
rx_sig = conv(qam, h_eq, 'same');
rx_sig = awgn(rx_sig, 15, 'measured');

figure();
plot(real(rx_sig), imag(rx_sig), '.', 'MarkerSize', 10);
grid on;
axis square;