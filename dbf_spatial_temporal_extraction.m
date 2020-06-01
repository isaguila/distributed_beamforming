close all;
clear;
clc;

rng(0, 'twister');

%--------------------------------------------------------------------------
%--------- Code based on paper by P. Sriploy and M. Uthansakul ------------
% Nonfeedback Distributed Beamforming using Spatial-Temporal Extraction
% Modified to include frequency offset and fading -------------------------
%--------------------------------------------------------------------------

%---------- System Parameters based on 802.11 set of protocols ------------
fc     = 5e9;       % Carrier Frequency
ppm    = 0.5;       % Carrier Frequency Offset in PPM
fs     = 20e6;      % Sampling rate
Ntx    = 4;         % Number of single-antenna transmitter nodes
snr_db = 10;        % Channel SNR in dB
e      = -1;        % Coefficient from paper, e = exp(1i*pi) = -1
A_NN   = ones(Ntx); % Matrix from paper

M      = 16;

if (Ntx == 4) % if 4 then use values from paper
    phase_deg = [-4.5; 108.1; -105.9; -25.5];
else          % Otherwise use random phase offset values
    phase_deg = -180 + 360*rand(Ntx, 1);
    phase_deg(1) = 0;
end
phase_rad = phase_deg*pi/180;

for index = 2:Ntx
    A_NN(index, index-1) = e;
end

mod_block = qam_modulator(M);
demod_block = qam_demodulator(M);

% Assume frequency offset ~U(-fmax, fmax)
fmax = 100; % Set to 0 for no frequency offset
df = -fmax + 2*fmax*rand(Ntx, 1);
df(1) = 0;

% Rayleigh flat fading (one tap channels)
h = 1/sqrt(2)*(randn(1, Ntx) + 1i*randn(1, Ntx)); h(1) = real(h(1));
% h = ones(1, Ntx); % Use all ones if no fading desired

mod_block = mod_block.get_syms(1000*log2(M));

tx = struct();
tx.bits = mod_block.bits;
tx.syms = mod_block.syms;
tx.sig = repmat(tx.syms, 1, Ntx)*diag(exp(1i*phase_rad));
sig_len = length(tx.sig);

%--------------- Apply Frequency Offset to Signal -------------------------
t = (0:sig_len-1).'/fs;
for ii = 1:Ntx
    tx.sig(:, ii) = tx.sig(:, ii).*exp(1i*2*pi*df(ii)*t);
end

%--------------- Apply Channel to Signal ----------------------------------
for ii = 1:Ntx
    tx.sig(:, ii) = conv(tx.sig(:, ii), h(ii), 'same');
end

%----------- Calculate N0 (noise density to use for AWGN) -----------------
p_sig_avg = mean(vecnorm(tx.sig))^2/sig_len;
p_noise_avg = p_sig_avg / 10^(snr_db/10);

%------------- Non-coherent Distributed Transmit Beamforming --------------
noise_samples = sqrt(p_noise_avg/2)*randn(sig_len, 2)*[1; 1i];
rx = struct();
rx.sig_noncoherent = sum(tx.sig, 2) + noise_samples;

%--------------------- Transmit signal N times ----------------------------
w_bb = sqrt(p_noise_avg/2)*(randn(sig_len, Ntx) + 1i*randn(sig_len, Ntx));
rx.sig = (A_NN*tx.sig.').' + w_bb;  % Received signal with noise

%------------------- Spatial-Temporal Extraction --------------------------
y = rx.sig*A_NN^(-1).';

%-------------------------- Phase Search ----------------------------------
phase_search_range = (-180:10:180)*pi/180;
phase_rad_hat = zeros(Ntx, 1);
for ii = 2:Ntx
    max_sig_pow = 0;
    for jj = 1:length(phase_search_range)
        phase_guess = phase_search_range(jj);
        sig_pow = norm(y(:, 1) + y(:, ii)*exp(1i*phase_guess))^2;
        if (sig_pow > max_sig_pow)
            max_sig_pow = sig_pow;
            phase_rad_hat(ii) = phase_guess;
        end
    end
end
phase_deg_hat = phase_rad_hat*180/pi;

%----------------- Precode signal with phase estimates --------------------
tx.sig_coherent = repmat(tx.syms, 1, Ntx)*diag(exp(1i*phase_rad_hat));

%----------------- Phase Offset (assume same as before) -------------------
tx.sig_coherent = tx.sig_coherent*diag(exp(1i*phase_rad));
sig_len = length(tx.sig_coherent);

%------- Apply Frequency Offset to Signal (assume same as before) ---------
t = (0:sig_len-1).'/fs;
for ii = 1:Ntx
    tx.sig_coherent(:, ii) = tx.sig_coherent(:, ii).*exp(1i*2*pi*df(ii)*t);
end

%----- Apply Channel to Signal (assume channel has not changed fading) ----
for ii = 1:Ntx
    tx.sig_coherent(:, ii) = conv(tx.sig_coherent(:, ii), h(ii), 'same');
end

%---------------- Coherent Distributed Transmit Beamforming ---------------
noise_samples = sqrt(p_noise_avg/2)*randn(sig_len, 2)*[1; 1i];
rx.sig_coherent = sum(tx.sig_coherent, 2) + noise_samples;

%--------------- Process the received signals -----------------------------
% Remove phase offset bias from LO and channel since all transmitters 
% align to the first antenna (i.e. constellation will be rotated if first
% antenna has phase offset other than 0)
rx.sig_coherent = rx.sig_coherent*exp(-1i*phase_rad(1));
rx.sig_coherent = rx.sig_coherent*exp(-1i*angle(h(1)));

%--------------- Equalize the signal --------------------------------------
% Need better parameters, this worsens performance!
% rx.sig_coherent = lms_equalizer(rx.sig_coherent, tx.syms, length(tx.syms), 3, 0.1);

demod_block = demod_block.demod_sig(rx.sig_coherent, sqrt(sig_len));
rx.bits = demod_block.bits;
rx.syms = demod_block.syms;

evm_db = 10*log10(1/length(rx.syms)*norm(rx.syms - tx.syms)^2);

phi_table = table(...
    phase_deg + wrapTo180(180/pi*reshape(angle(h), size(phase_deg))), ...
    phase_deg_hat, ...
    'VariableNames', {'Total Phase Offset (degrees)', 'Phase Compensation (degrees)'});
disp(phi_table);

figure();
plot(real(tx.syms), imag(tx.syms), 'bx', 'MarkerSize', 25, 'LineWidth', 2);
hold on;
plot(real(rx.syms), imag(rx.syms), 'r.', 'MarkerSize', 10);
legend('Tx', 'Rx');
xlabel('In-Phase');
ylabel('Quadrature');
title('Ideal Constellation and Received Constellation');
grid on;
xlim([-2 2]);
ylim([-2 2]);
axis square;

bit_errs = sum(bitxor(tx.bits, rx.bits));
fprintf('Bit Errors:        %d\n', bit_errs);
fprintf('Channel SNR:       %.2f dB\n', snr_db);
fprintf('Measured SNR:      %.2f dB\n', -evm_db);
fprintf('Theoretical Gain:  %.2f dB\n', 20*log10(Ntx));
fprintf('Measured Gain:     %.2f dB\n', -evm_db - snr_db);
