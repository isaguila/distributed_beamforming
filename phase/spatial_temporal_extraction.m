close all;
clear;
clc;

rng(0, 'twister');

%--------------------------------------------------------------------------
%--------- Code based on paper by P. Sriploy and M. Uthansakul ------------
% Nonfeedback Distributed Beamforming using Spatial-Temporal Extraction
%--------------------------------------------------------------------------

%---------- System Parameters based on 802.11 set of protocols ------------
fs     = 20e6;      % Sampling rate
Ntx    = 20;         % Number of single-antenna transmitter nodes
snr_db = 15;        % Channel SNR in dB
e      = -1;        % Coefficient from paper, e = exp(1i*pi) = -1
A_NN   = ones(Ntx); % Matrix from paper

sig_model_list = {...
    'Tone'; ...          % 1
    'BPSK'; ...          % 2
    '4-QAM'; ...         % 3
    '16-QAM'; ...        % 4
    '64-QAM'; ...        % 5
    '256-QAM'; ...       % 6
    '1024-QAM'; ...      % 7
    'OFDM BPSK'; ...     % 8
    'OFDM 4-QAM'; ...    % 9
    'OFDM 16-QAM'; ...   % 10
    'OFDM 64-QAM'; ...   % 11
    'OFDM 256-QAM'; ...  % 12
    'OFDM 1024-QAM'; ... % 13
    };

sig_model = sig_model_list{7};

if (Ntx == 4) % if 4 then use values from paper
    phase_deg = [-4.5; 108.1; -105.9; -25.5];
else          % Otherwise use random phase offset values
    phase_deg = -180 + 360*rand(Ntx, 1);
end
phase_rad = phase_deg*pi/180;

for index = 2:Ntx
    A_NN(index, index-1) = e;
end

tx = struct();
rx = struct();

tx_bl = struct();
rx_bl = struct();


switch sig_model
    case 'Tone'
        fc = fs/32;
        tone_len = 10000;
        t = (0:tone_len-1)'/fs;
        tx.sig = complex(zeros(tone_len, Ntx));
        for index = 1:Ntx
            tx.sig(:, index) = exp(1i*(2*pi*fc*t + phase_rad(index)));
        end
    case 'BPSK'
        num_psk_syms = 10000;
        tx.bits = randi([0 1], num_psk_syms, 1);
        tx.syms = 2*tx.bits-1;
        tx.sig = repmat(tx.syms, 1, Ntx)*diag(exp(1i*phase_rad));
    case '4-QAM'
        M = 4;
        num_qam_syms = 1000*log2(M);
        tx.bits = randi([0 1], log2(M)*num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.sig = repmat(tx.syms, 1, Ntx)*diag(exp(1i*phase_rad));
    case '16-QAM'
        M = 16;
        num_qam_syms = 1000*log2(M);
        tx.bits = randi([0 1], log2(M)*num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.sig = repmat(tx.syms, 1, Ntx)*diag(exp(1i*phase_rad));
    case '64-QAM'
        M = 64;
        num_qam_syms = 1000*log2(M);
        tx.bits = randi([0 1], log2(M)*num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.sig = repmat(tx.syms, 1, Ntx)*diag(exp(1i*phase_rad));
    case '256-QAM'
        M = 256;
        num_qam_syms = 1000*log2(M);
        tx.bits = randi([0 1], log2(M)*num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.sig = repmat(tx.syms, 1, Ntx)*diag(exp(1i*phase_rad));
    case '1024-QAM'
        M = 1024;
        num_qam_syms = 1000*log2(M);
        tx.bits = randi([0 1], log2(M)*num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.sig = repmat(tx.syms, 1, Ntx)*diag(exp(1i*phase_rad));
    case 'OFDM BPSK'
        fft_len = 64;
        guard_len = 16;
        ofdm_sym_len = fft_len + guard_len;
        num_psk_syms = 1000*fft_len;
        tx.bits = randi([0 1], num_psk_syms, 1);
        tx.syms = 2*tx.bits-1;
        tx.syms_matrix = reshape(tx.syms, fft_len, length(tx.syms)/fft_len);
        tx.syms_matrix = ifft(ifftshift(tx.syms_matrix), fft_len);
        tx.guard_int = tx.syms_matrix(end-guard_len+1:end, :);
        tx.ofdm_syms = [tx.guard_int; tx.syms_matrix];
        tx.sig = repmat(reshape(tx.ofdm_syms, [], 1), 1, Ntx)*diag(exp(1i*phase_rad));
    case 'OFDM 4-QAM'
        M = 4;
        fft_len = 64;
        guard_len = 16;
        ofdm_sym_len = fft_len + guard_len;
        num_qam_syms = 1000*log2(M)*fft_len;
        tx.bits = randi([0 1], num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.syms_matrix = reshape(tx.syms, fft_len, length(tx.syms)/fft_len);
        tx.syms_matrix = ifft(ifftshift(tx.syms_matrix), fft_len);
        tx.guard_int = tx.syms_matrix(end-guard_len+1:end, :);
        tx.ofdm_syms = [tx.guard_int; tx.syms_matrix];
        tx.sig = repmat(reshape(tx.ofdm_syms, [], 1), 1, Ntx)*diag(exp(1i*phase_rad));
    case 'OFDM 16-QAM'
        M = 16;
        fft_len = 64;
        guard_len = 16;
        ofdm_sym_len = fft_len + guard_len;
        num_qam_syms = 1000*log2(M)*fft_len;
        tx.bits = randi([0 1], num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.syms_matrix = reshape(tx.syms, fft_len, length(tx.syms)/fft_len);
        tx.syms_matrix = ifft(ifftshift(tx.syms_matrix), fft_len);
        tx.guard_int = tx.syms_matrix(end-guard_len+1:end, :);
        tx.ofdm_syms = [tx.guard_int; tx.syms_matrix];
        tx.sig = repmat(reshape(tx.ofdm_syms, [], 1), 1, Ntx)*diag(exp(1i*phase_rad));
    case 'OFDM 64-QAM'
        M = 64;
        fft_len = 64;
        guard_len = 16;
        ofdm_sym_len = fft_len + guard_len;
        num_qam_syms = 1000*log2(M)*fft_len;
        tx.bits = randi([0 1], num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.syms_matrix = reshape(tx.syms, fft_len, length(tx.syms)/fft_len);
        tx.syms_matrix = ifft(ifftshift(tx.syms_matrix), fft_len);
        tx.guard_int = tx.syms_matrix(end-guard_len+1:end, :);
        tx.ofdm_syms = [tx.guard_int; tx.syms_matrix];
        tx.sig = repmat(reshape(tx.ofdm_syms, [], 1), 1, Ntx)*diag(exp(1i*phase_rad));
    case 'OFDM 256-QAM'
        M = 256;
        fft_len = 64;
        guard_len = 16;
        ofdm_sym_len = fft_len + guard_len;
        num_qam_syms = 1000*log2(M)*fft_len;
        tx.bits = randi([0 1], num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.syms_matrix = reshape(tx.syms, fft_len, length(tx.syms)/fft_len);
        tx.syms_matrix = ifft(ifftshift(tx.syms_matrix), fft_len);
        tx.guard_int = tx.syms_matrix(end-guard_len+1:end, :);
        tx.ofdm_syms = [tx.guard_int; tx.syms_matrix];
        tx.sig = repmat(reshape(tx.ofdm_syms, [], 1), 1, Ntx)*diag(exp(1i*phase_rad));
    case 'OFDM 1024-QAM'
        M = 1024;
        fft_len = 64;
        guard_len = 16;
        ofdm_sym_len = fft_len + guard_len;
        num_qam_syms = 1000*log2(M)*fft_len;
        tx.bits = randi([0 1], num_qam_syms, 1);
        tx.syms = qammod(tx.bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        tx.syms_matrix = reshape(tx.syms, fft_len, length(tx.syms)/fft_len);
        tx.syms_matrix = ifft(ifftshift(tx.syms_matrix), fft_len);
        tx.guard_int = tx.syms_matrix(end-guard_len+1:end, :);
        tx.ofdm_syms = [tx.guard_int; tx.syms_matrix];
        tx.sig = repmat(reshape(tx.ofdm_syms, [], 1), 1, Ntx)*diag(exp(1i*phase_rad));
    otherwise
        error('Not implemented yet!');
end

tx_bl.sig = tx.sig(:, 1);

%----------- Calculate N0 (noise density to use for AWGN) -----------------
sig_len = length(tx.sig);
p_sig_avg = mean(vecnorm(tx.sig))^2/sig_len;
p_noise_avg = p_sig_avg / 10^(snr_db/10);

%---------------- baseline

%------------- Non-coherent Distributed Transmit Beamforming --------------
noise_samples = sqrt(p_noise_avg/2)*randn(sig_len, 2)*[1; 1i];
rx.sig_noncoherent = sum(tx.sig, 2) + noise_samples;
rx_bl.sig = tx_bl.sig + noise_samples;

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
tx.sig_precoded = tx.sig*diag(exp(1i*phase_rad_hat));

%---------------- Coherent Distributed Transmit Beamforming ---------------
noise_samples = sqrt(p_noise_avg/2)*randn(sig_len, 2)*[1; 1i];
rx.sig_coherent = (sum(tx.sig_precoded, 2) + noise_samples);

%--------------- Process the received signals -----------------------------
switch sig_model
    case 'Tone'
        % Do nothing
    case 'BPSK'
        rx.syms = rx.sig_coherent*exp(-1i*phase_rad(1));
        rx.syms = rx.syms*sqrt(sig_len)/norm(rx.syms);
        rx.bits = nan(size(tx.bits));
        rx.bits(rx.syms < 0) = 0;
        rx.bits(rx.syms > 0) = 1;
    case {'4-QAM', '16-QAM', '64-QAM', '256-QAM', '1024-QAM'}
        rx_bl.syms = rx_bl.sig*exp(-1i*phase_rad(1));
        rx_bl.syms = rx_bl.syms*sqrt(sig_len)/norm(rx_bl.syms);
        rx_bl.bits = qamdemod(rx_bl.syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);
        
        rx.syms = rx.sig_coherent*exp(-1i*phase_rad(1));
        rx.syms = rx.syms*sqrt(sig_len)/norm(rx.syms);
        rx.bits = qamdemod(rx.syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    case 'OFDM BPSK'
        rx.ofdm_syms = reshape(rx.sig_coherent, ofdm_sym_len, sig_len/ofdm_sym_len);
        rx.ofdm_syms = rx.ofdm_syms*exp(-1i*phase_rad(1));
        rx.syms_matrix = rx.ofdm_syms(guard_len+1:end, :);
        rx.syms_matrix  = fftshift(fft(rx.syms_matrix, fft_len));
        rx.syms = reshape(rx.syms_matrix, [], 1);
        rx.syms = rx.syms*sqrt(sig_len)/norm(rx.syms);
        rx.bits = nan(size(tx.bits));
        rx.bits(rx.syms < 0) = 0;
        rx.bits(rx.syms > 0) = 1;
    case {'OFDM 4-QAM', 'OFDM 16-QAM', 'OFDM 64-QAM', 'OFDM 256-QAM', 'OFDM 1024-QAM'}
        rx.ofdm_syms = reshape(rx.sig_coherent, ofdm_sym_len, sig_len/ofdm_sym_len);
        rx.ofdm_syms = rx.ofdm_syms*exp(-1i*phase_rad(1));
        rx.syms_matrix = rx.ofdm_syms(guard_len+1:end, :);
        rx.syms_matrix  = fftshift(fft(rx.syms_matrix, fft_len));
        rx.syms = reshape(rx.syms_matrix, [], 1);
        rx.syms = rx.syms*norm(tx.syms)/norm(rx.syms);
        rx.bits = qamdemod(rx.syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    otherwise
        error('Not implemented yet!');
end

evm_db = 10*log10(1/length(rx.syms)*norm(rx.syms - tx.syms)^2);
evm_siso = 10*log10(1/length(rx.syms)*norm(rx_bl.syms - tx.syms)^2);

fsize = 14;

figure();
plot(real(rx.syms), imag(rx.syms), 'r.', 'MarkerSize', 10);
xlabel('In-Phase', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Quadrature', 'FontSize', fsize, 'Interpreter', 'Latex');
title('Spatial Temporal Extraction', 'FontSize', fsize, 'Interpreter', 'Latex');
legend(['EVM = ', num2str(evm_db), ' dB'], 'FontSize', fsize, 'Interpreter', 'Latex');
grid on;
xlim([-2 2]);
ylim([-2 2]);
axis square;
set(gcf, 'color', 'w');

fn = sprintf('/Users/ivan/Documents/Thesis/figures/N%d_ste_perf_%dqam', Ntx, M);
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');


figure();
plot(real(rx_bl.syms), imag(rx_bl.syms), 'r.', 'MarkerSize', 10);
xlabel('In-Phase', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Quadrature', 'FontSize', fsize, 'Interpreter', 'Latex');
title('Baseline - SISO', 'FontSize', fsize, 'Interpreter', 'Latex');
legend(['EVM = ', num2str(evm_siso), ' dB'], 'FontSize', fsize, 'Interpreter', 'Latex');
grid on;
xlim([-2 2]);
ylim([-2 2]);
axis square;
set(gcf, 'color', 'w');

fn = sprintf('/Users/ivan/Documents/Thesis/figures/N%d_ste_siso_perf_%dqam', Ntx, M);
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');


bit_errs = sum(bitxor(tx.bits, rx.bits));
fprintf('Bit Errors:       %d\n', bit_errs);
fprintf('Channel SNR:      %.2f dB\n', snr_db);
fprintf('Measured SNR:     %.2f dB\n', -evm_db);
fprintf('Theoretical Gain: %.2f dB\n', 20*log10(Ntx));
fprintf('Measured Gain:    %.2f dB\n', -evm_db - snr_db);