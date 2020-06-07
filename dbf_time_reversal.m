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

Ntx = 4;
M = 4;
snr_db = 20;

mod_block = qam_modulator(M);
demod_block = qam_demodulator(M);
mod_block = mod_block.get_syms(5000*log2(M));
scale_factor = vecnorm(mod_block.syms).';

tx_tr = struct();
rx_tr = struct();

tx_tr.sig    = repmat(mod_block.syms, 1, Ntx);

L = 4000;
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
end
 

sig_len = length(tx_tr.sig);
p_sig = mean(vecnorm(tx_tr.sig).^2);
p_noise = p_sig / 10^(snr_db);
noise = sqrt(p_noise/2)*randn(sig_len, 2)*[1; 1i];
rx_tr.sig = sum(tx_tr.sig, 2) + noise;


demod_block = demod_block.demod_sig(rx_tr.sig, scale_factor);
rx_tr.syms = reshape(demod_block.syms, [], 1);
figure();
plot(real(rx_tr.syms), imag(rx_tr.syms), '.', 'MarkerSize', 10);
title('Time Reversal');
grid on;
xlim([-2, 2]);
ylim([-2, 2]);
axis square;