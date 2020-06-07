close all;
clear;
clc;

rng(0);

%--------------------------------------------------------------------------
%------------ Distributed Beamforming with Time Reversal ------------------
%------------------------- Code based on papers ---------------------------
% Multipath effects on time reversal OFDM communications between wireless sensors
% by Z. Chen, Y. Zhao, D. Zhao
% 
% Time Reversal Techniques for Wireless Communications 
% by P. Kyritsi, G. Papanicolaou, P. Eggers, A. Oprea
%--------------------------------------------------------------------------
Ntx = 10;
M = 4;
snr_db = 15;
L = 1e4;


mod_block = qam_modulator(M);
demod_block = qam_demodulator(M);
mod_block = mod_block.get_syms(5000*log2(M));
scale_factor = vecnorm(mod_block.syms).';

tx_tr = struct();
rx_tr = struct();

tx_bl = struct();
rx_bl = struct();

tx_ntr = struct();
rx_ntr = struct();

tx_tr.sig  = repmat(mod_block.syms, 1, Ntx);
tx_ntr.sig = repmat(mod_block.syms, 1, Ntx);
tx_bl.sig  = mod_block.syms;

h = complex(zeros(L, Ntx));
h_tr = complex(zeros(L, Ntx));

for ii = 1:Ntx
    % Create the channel and time reversed version
    h(:, ii) = randn(L, 2)*[1; 1i];
    h_tr(:, ii) = flip(conj(h(:, ii)));
    
    % Precompensate TR 
    tx_tr.sig(:, ii) = conv(tx_tr.sig(:, ii), h_tr(:, ii), 'same');
     
    % Apply channel
    tx_tr.sig(:, ii) = conv(tx_tr.sig(:, ii), h(:, ii), 'same');
    tx_ntr.sig(:, ii) = conv(tx_ntr.sig(:, ii), h(:, ii), 'same');
end

% Choose a channel for baseline SISO and apply it
h_bl = randn(L, 2)*[1; 1i];
tx_bl.sig = conv(tx_bl.sig, h_bl, 'same');

 
sig_len = length(tx_tr.sig);
p_sig = mean(vecnorm(tx_tr.sig).^2)/sig_len;
p_noise = p_sig / 10^(snr_db/10);
noise = sqrt(p_noise/2)*randn(sig_len, 2)*[1; 1i];
rx_tr.sig = sum(tx_tr.sig, 2) + noise;

sig_len = length(tx_ntr.sig);
p_sig = mean(vecnorm(tx_ntr.sig).^2)/sig_len;
p_noise = p_sig / 10^(snr_db/10);
noise = sqrt(p_noise/2)*randn(sig_len, 2)*[1; 1i];
rx_ntr.sig = sum(tx_ntr.sig, 2) + noise;

sig_len = length(tx_bl.sig);
p_sig = mean(vecnorm(tx_bl.sig).^2)/sig_len;
p_noise = p_sig / 10^(snr_db/10);
noise = sqrt(p_noise/2)*randn(sig_len, 2)*[1; 1i];
rx_bl.sig = sum(tx_bl.sig, 2) + noise;

% phi = angle(rx_tr.sig.*conj(mod_block.syms)
% phi = angle(sum(rx_tr.sig.*conj(mod_block.syms)));
% rx_tr.sig = rx_tr.sig.*exp(-1i*phi);
demod_block = demod_block.demod_sig(rx_tr.sig, scale_factor);
rx_tr.syms = reshape(demod_block.syms, [], 1);
rx_tr.evm = get_evm_qam(rx_tr.syms, M);
rx_tr.bits = demod_block.bits;
rx_tr.ber = sum(bitxor(rx_tr.bits, mod_block.bits));

demod_block = demod_block.demod_sig(rx_ntr.sig, scale_factor);
rx_ntr.syms = reshape(demod_block.syms, [], 1);
rx_ntr.evm = get_evm_qam(rx_ntr.syms, M);
rx_ntr.bits = demod_block.bits;
rx_ntr.ber = sum(bitxor(rx_ntr.bits, mod_block.bits));

demod_block = demod_block.demod_sig(rx_bl.sig, scale_factor);
rx_bl.syms = reshape(demod_block.syms, [], 1);
rx_bl.evm = get_evm_qam(rx_bl.syms, M);
rx_bl.bits = demod_block.bits;
rx_bl.ber = sum(bitxor(rx_bl.bits, mod_block.bits));

figure();
plot(real(rx_tr.syms), imag(rx_tr.syms), '.', 'MarkerSize', 10);
title('Time Reversal');
grid on;
xlim([-2, 2]);
ylim([-2, 2]);
axis square;

figure();
plot(real(rx_ntr.syms), imag(rx_ntr.syms), '.', 'MarkerSize', 10);
title('Non Time Reversal');
grid on;
xlim([-2, 2]);
ylim([-2, 2]);
axis square;

figure();
plot(real(rx_bl.syms), imag(rx_bl.syms), '.', 'MarkerSize', 10);
title('Baseline');
grid on;
xlim([-2, 2]);
ylim([-2, 2]);

disp(rx_tr.evm);
disp(rx_ntr.evm);
disp(rx_bl.evm);

% scatterplot(awgn(qammod(randi([0 M-1], 1000, 1), M, 'UnitAveragePower', true), -rx_tr.evm, 'measured'))
