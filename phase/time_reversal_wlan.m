close all;
clear;
clc;

rng(0);

tgn = wlanTGnChannel(...
    'PathGainsOutputPort', true, ...
    'DelayProfile', 'Model-E');

bits = randi([0 1], 10000, 1);
x = qammod(bits, 4, 'InputType', 'bit', 'UnitAveragePower', 1);
snr = 20;

[y, path_gains] = tgn(x);

h = path_gains(1, :).';
phase_response = angle(h);
mag_response = 20*log10(abs(h));

h = [1; 0; 0; 0; 0.5i];

x_tr = conv(x, conj(flip(h)));
% [y_tr, chan_coefs] = tgn(x_tr);
y_tr = conv(x_tr, h);
y_tr = y_tr*norm(x)/norm(y_tr);

p_noise = 1/10^(snr/10);
noise1 = sqrt(p_noise/2)*randn(length(y), 2)*[1; 1i];
noise2 = sqrt(p_noise/2)*randn(length(y_tr), 2)*[1; 1i];
y = y + noise1;
y_tr = y_tr + noise2;

figure(1);
plot(real(y), imag(y), '.', 'MarkerSize', 10, 'DisplayName', 'Before Equalization');
grid on;
xlabel('In-Phase');
ylabel('Quadrature');
title('y(t)');

figure(2);
plot(real(y_tr), imag(y_tr), '.', 'MarkerSize', 10, 'DisplayName', 'Before Equalization');
grid on;
xlabel('In-Phase');
ylabel('Quadrature');
title('y_{tr}(t)');

numPoints = length(y);
numTaps = 20;		% channel order
Mu = 0.01;          % iteration step size

w1 = [];
z1 = [];
in1 = [];
e1 = []; % error, final result to be computed
w1 = zeros(numTaps+1,1) + 1i*zeros(numTaps+1,1);

% LMS Adaptation
for n  = numTaps+1 : numPoints
    % select part of training input
    in1 = x(n : -1 : n-numTaps) ;
    z1(n) = w1'*in1;
    % compute error
    e1(n) = y(n)-z1(n);
    % update taps
    w1 = w1 + Mu*( real(e1(n)*conj(in1)) - 1i*imag(e1(n)*conj(in1)) );
end

numPoints = length(x);
numTaps = 9;		% channel order
Mu = 0.1;           % iteration step size

w2 = [];
zz = [];
in2 = [];
e2 = []; % error, final result to be computed
w2 = zeros(numTaps+1,1) + 1i*zeros(numTaps+1,1);

% LMS Adaptation
for n  = numTaps+1 : numPoints
    % select part of training input
    in2 = x(n : -1 : n-numTaps) ;
    z2(n) = w2'*in2;
    % compute error
    e2(n) = y_tr(n)-z2(n);
    % update taps
    w2 = w2 + Mu*( real(e2(n)*conj(in2)) - 1i*imag(e2(n)*conj(in2)) );
end

figure(1);
hold on;
plot(real(z1), imag(z1), '.', 'MarkerSize', 10, 'DisplayName', 'After Equalization');
legend('show');
axis square;

figure(2);
hold on;
plot(real(z2), imag(z2), '.', 'MarkerSize', 10, 'DisplayName', 'After Equalization');
legend('show');
axis square;


