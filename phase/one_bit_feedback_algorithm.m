close all;
clear;
clc;

fsize = 14;
fig_ctr = 0;
Ntx = 5;

fs = 10e3;
fc = 100;

t = 0:1/fs:10/fc; t = t(:);
t_ms = t*1e3;
M = length(t);
phi = -pi + 2*pi*rand(Ntx, 1);

x = zeros(M, 1);
figure();
subplot(2, 2, 1);
for ii = 1:Ntx
    x(:, ii) = exp(1i*(2*pi*fc*t + phi(ii)));
    plot(t_ms, real(x(:, ii)), 'DisplayName', ['$\phi_', num2str(ii), '$ = ', num2str(phi(ii)*180/pi), '{$^\circ$}']);
    xlabel('Time (ms)', 'FontSize', fsize, 'Interpreter', 'Latex');
    ylabel('Amplitude', 'FontSize', fsize, 'Interpreter', 'Latex');
    hold on;
    grid on;
    title('Starting Signals', 'FontSize', fsize, 'Interpreter', 'Latex');
end
set(gcf, 'color', 'w');
legend('show', 'FontSize', 12, 'Interpreter', 'Latex');
hold off;

iter = 100*Ntx;
mu = 0.1;
phi_hat = zeros(iter, Ntx);
phi_hat(1, :) = phi;

energy_vec = zeros(iter, 1);
energy_vec(1) = norm(sum(x, 2))^2/M;

sig_energy = energy_vec;

for idx = 2:iter
    rand_weight = mu*2*pi*rand(Ntx, 1);
    
    y = x*diag(exp(1i*rand_weight));
    
    energy_vec(idx) = norm(sum(y, 2))^2/M; 
    
    if energy_vec(idx) > energy_vec(idx-1)
        x = y;
        phi_hat(idx, :) = rand_weight;
    else
        phi_hat(idx, :) = phi_hat(idx-1, :);
    end
    
    sig_energy(idx) = norm(sum(x, 2))^2/M; 
    z = 1;
end

fig_ctr = fig_ctr + 1;
figure(fig_ctr);
subplot(2, 2, 2);
for ii = 1:Ntx
    plot(t_ms, real(x(:, ii)), 'DisplayName', ['$\phi_', num2str(ii), '$ = ', num2str(phi_hat(end, ii)*180/pi), '{$^\circ$}']);
    hold on;
    grid on;
    title('Ending Signals', 'FontSize', fsize, 'Interpreter', 'Latex');
    xlabel('Time (ms)', 'FontSize', fsize, 'Interpreter', 'Latex');
end
set(gcf, 'color', 'w');
legend('show', 'FontSize', fsize, 'Interpreter', 'Latex');
hold off;

subplot(2, 2, 3:4);
plot(sig_energy/Ntx^2);
grid on;
title('Signal Energy', 'FontSize', fsize, 'Interpreter', 'Latex');
xlabel('Iteration Number', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Normalized Energy', 'FontSize', fsize, 'Interpreter', 'Latex');
set(gcf, 'color', 'w');

fn = '/Users/ivan/Documents/Thesis/figures/one_bit_feedback_energy';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');

final_energy = norm(sum(x, 2))^2/M;
fprintf('Final Energy = %.3f\n', final_energy);
fprintf('Percent achieved = %.3f\n', final_energy/Ntx^2);

fig_ctr = fig_ctr + 1;
figure(fig_ctr);
for ii = 1:Ntx
   polarplot(repmat(phi_hat(1, ii), 2, 1), [1; 0], 'LineWidth', 2, 'DisplayName', ['$\phi_', num2str(ii), ' = $', num2str(phi_hat(end, ii)*180/pi), '{$^\circ$}']);
   hold on;
end
legend('show', 'FontSize', 12, 'Interpreter', 'Latex');
set(gcf, 'color', 'w');
fn = '/Users/ivan/Documents/Thesis/figures/one_bit_feedback_polar_start';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');

fig_ctr = fig_ctr + 1;
figure(fig_ctr);
for ii = 1:Ntx
   polarplot(repmat(phi_hat(end, ii), 2, 1), [1; 0], 'LineWidth', 2, 'DisplayName', ['$\phi_', num2str(ii), ' = $', num2str(phi_hat(end, ii)*180/pi), '{$^\circ$}']);
   hold on;
end
legend('show', 'FontSize', 12, 'Interpreter', 'Latex');
set(gcf, 'color', 'w');
fn = '/Users/ivan/Documents/Thesis/figures/one_bit_feedback_polar_end';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');

%%
figure();
rho_ = repmat(phi_hat(end, :), 2, 1);
polarplot(rho_, [zeros(1, Ntx); ones(1, Ntx)], 'LineWidth', 2);
legend({'1', '2', '3', '4', '5'}, 'FontSize', 12, 'Interpreter', 'Latex');