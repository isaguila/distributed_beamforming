close all;
clear;
clc;

rng(0, 'twister');

%--------------------------------------------------------------------------
%------- Code based on paper by D. R. Brown III, Y. Liao, N. Fox ----------
% Low-Complexity Real-Time Single Tone Phase and Frequency Estimation
%--------------------------------------------------------------------------

fs = 16e3;
N = 1000;
t_vec = (0:N-1).'/fs;

alpha = 0.3;

fc = 1010;
fc = fc + (1030 - 1010)*rand(1, 1);
phi_0 = 0;
phi_0 = phi_0 + -pi + 2*pi*rand(1, 1);
x = exp(1i*(2*pi*fc*t_vec + phi_0));
x = awgn(x, 10, 'measured');

zc_detect = false;

l = [-1, 0, 1];
l_p = l;
l_n = l;

l_len = numel(l);

k_p = zeros(1, l_len);
k_n = zeros(1, l_len);

index = 1;

% Neg Crossing occur at
% 1 -> 2
% 3 -> 2
% Pos Crossing occur at
% 2 -> 1
% 4 -> 1

% Real part of signal
% + crossing must occur at phase k*2*pi - pi/2
% - crossing must occur at phase k*2*pi + pi/2

% Imag part of signal
% + crossing must occur at phase l*2*pi
% - crossing must occur at phase l*2*pi + pi

%-- Init state variables
A = 0;
B = 0;
C = 0;
D = zeros(1, l_len);
E = zeros(1, l_len);

t_array = zeros(N, 1);
phase_array = zeros(N, l_len);

phi_hat = zeros(N, l_len);
omega_hat = zeros(N, l_len);

y = [real(x), imag(x)];
state = zeros(1, l_len);
temp = zeros(1, l_len);
current_state = 0;

figure();
plot(alpha*ones(N, 1), 'b', 'LineWidth', 2);
hold on;
plot(-alpha*ones(N, 1), 'b', 'LineWidth', 2);
plot(y(:, 1), 'b');
plot(y(:, 2), 'r');
legend('','','Real', 'Imag');
grid on;
hold off;

for ii = 2:N
    for iq = 1:2
        val = y(ii, iq);
        zc_detect = false;
        current_state = state(iq);
        switch current_state
            case 1
                if abs(val) <= alpha
                    current_state = 3;
                    temp(iq) = ii;
                elseif val < -alpha
                    % --------- Negative Zero Crossing (1 -> 2) -----------
                    current_state = 2;
                    zc_detect = true;
                    
                    if iq == 1
                        phase_i = k_n.*2*pi + pi/2;
                        k_n = k_n + ones(1, l_len);
                    else
                        phase_i = l_n.*2*pi + pi;
                        l_n = l_n + ones(1, l_len);
                    end
                end
            case 2
                if abs(val) <= alpha
                    current_state = 4;
                    temp(iq) = ii;
                elseif val > alpha
                    % -------- Positive Zero Crossing (2 -> 1) ------------
                    current_state = 1;
                    zc_detect = true;
                    
                    if iq == 1
                        phase_i = k_p.*2*pi - pi/2;
                        k_p = k_p + ones(1, l_len);
                    else
                        phase_i = l_p.*2*pi;
                        l_p = l_p + ones(1, l_len);
                    end
                end
            case 3
                if val > alpha
                    current_state = 1;
                    temp(iq) = ii;
                elseif val < -alpha
                    % ---------- Negative Zero Crossing (3 -> 2) ----------
                    current_state = 2;
                    zc_detect = true;
                    
                    if iq == 1
                        phase_i = k_n.*2*pi + pi/2;
                        k_n = k_n + ones(1, l_len);
                    else
                        phase_i = l_n.*2*pi + pi;
                        l_n = l_n + ones(1, l_len);
                    end
                end
            case 4
                if val < -alpha
                    current_state = 2;
                    temp(iq) = ii;
                elseif val > alpha
                    %-------- Positive Zero Crossing (4 -> 1) -------------
                    current_state = 1;
                    zc_detect = true;
                    
                    if iq == 1
                        phase_i = k_p.*2*pi - pi/2;
                        k_p = k_p + ones(1, l_len);
                    else
                        phase_i = l_p.*2*pi;
                        l_p = l_p + ones(1, l_len);
                    end
                end
            otherwise
                current_state = 0;
                temp(iq) = ii;
                if val > alpha
                    current_state = 1;
                elseif val < -alpha
                    current_state = 2;
                end
        end
        state(iq) = current_state;
        
        if zc_detect == true
            % y0 = mt0 + b
            % yi = mti + b
            % y0 - yi = mt0 - mti = -yi
            % -yi/m = t0 - ti
            % t0 = ti - yi/m
            
            t_curr_state = ii;
            t_prev_state = temp(iq)-1;
            
            %-- Update variables
            m = (val - y(t_prev_state, iq))/(t_curr_state - t_prev_state);
            t_i = t_curr_state - val/m;
            t_i = t_i/fs;
            
            t_array(index) = t_i;
            phase_array(index, :) = phase_i;
            
            %-- Update State Variables
            A = A + t_i^2;
            B = B + 1;
            C = C + t_i;
            D = D - t_i*phase_i;
            E = E - phase_i;
            phi_hat(index, :) = -(E.*A - C.*D)./(B*A - C^2);
            omega_hat(index, :) = -C*phi_hat(index, :)./A - D./A;
            
            index = index + 1;
        end
    end
    
end

index = index - 1;
phase_array(index+1:end, :) = [];
t_array(index+1:end) = [];
phi_hat(index+1:end, :) = [];
omega_hat(index+1:end, :) = [];

omega_est = omega_hat(index, :);
phi_est = phi_hat(index, :);

epsilon = zeros(1, l_len);
mse = zeros(1, l_len);

figure();
for ii = 1:l_len
    epsilon(ii) = norm(exp(1i*(omega_est(ii)*t_vec + phi_est(ii))) - x)^2;
    mse(ii) = fc - omega_est(ii)/2/pi;
    
    fprintf('Fc = %.2f Hz, Phi = %.2f\n', fc, phi_0*180/pi);
    fprintf('Fe = %.2f Hz, Phi = %.2f\n', omega_est(ii)/2/pi, phi_est(ii)*180/pi);
    fprintf('Epsilon = %.2f\n', epsilon(ii));
    fprintf('Mean Square Error = %.2f\n', mse(ii));
    fprintf('\n');
    
    title_str = sprintf('l = %d, Epsilon = %.4f, RMSE = %.4f', l(ii), epsilon(ii), mse(ii));
    subplot(l_len, 1, ii);
    plot(t_vec*1e3, real(x), 'r');
    hold on;
    plot(t_vec*1e3, real(exp(1i*(omega_est(ii)*t_vec + phi_est(ii)))), 'b');
    legend('Actual', 'Estimated');
    title(title_str);
    grid on;
    xlabel('time (ms)');
    hold off;
end
