close all;
clear;
clc;

%--------------------------------------------------------------------------
%------- Code based on paper by R. Mudumbai, G. Barriac, U. Madhow --------
% On the Feasibility of Distributed Beamforming in Wireless Networks
%--------------------------------------------------------------------------
fsize = 14;

Ntx       = 2;          % Number of transmitter nodes
c         = 299792458;  % speed of light in m/s
% phi_ideal = wrapToPi(2*pi*d*fc/c)*180/pi; % ideal phase in degrees

%--------------------------------------------------------------------------
%------------------------ Frequency Errors --------------------------------
%---------------------------- Fc = 915 MHz --------------------------------
fc = 915e6;         % carrier frequency in Hz
d  = 500;           % Distance to receiver
df = -1e4:500:1e4;  % Relative frequency offset
phi_freq_error = wrapToPi(2*pi*d.*(fc + df)/c)*180/pi;

figure();
plot(1e-3*df, phi_freq_error, 'LineWidth', 2, 'DisplayName', ['d = ', num2str(d), ' m']);
hold on;
title(['Phase vs Frequency Error at $F_c$ = ', num2str(fc/1e6), ' MHz'], 'FontSize', fsize, 'Interpreter', 'Latex');
xlabel('Frequency Offset (kHz)', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Phase (degrees)', 'FontSize',fsize, 'Interpreter', 'Latex');
grid on;

d  = 1e3;         % Distance to receiver
phi_freq_error = wrapToPi(2*pi*d.*(fc + df)/c)*180/pi;
plot(1e-3*df, phi_freq_error, 'LineWidth', 2, 'DisplayName', ['d = ', num2str(d), ' m']);

d  = 1e4;         % Distance to receiver
phi_freq_error = wrapToPi(2*pi*d.*(fc + df)/c)*180/pi;
plot(1e-3*df, phi_freq_error, 'LineWidth', 2, 'DisplayName', ['d = ', num2str(d), ' m']);
set(gcf, 'color', 'w');

legend('show', 'FontSize', fsize);

fn = '/Users/Ivan/Documents/Thesis/figures/propagation_phase_vs_freq_915MHz';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');

%------------------------------ Fc = 5 GHz --------------------------------

fc = 5e9;           % carrier frequency in Hz
d  = 500;           % Distance to receiver
df = -1e4:500:1e4;  % Relative frequency offset
phi_freq_error = wrapToPi(2*pi*d.*(fc + df)/c)*180/pi;

figure();
plot(1e-3*df, phi_freq_error, 'LineWidth', 2, 'DisplayName', ['d = ', num2str(d), ' m']);
hold on;
title(['Phase vs Frequency Error at $F_c$ = ', num2str(fc/1e6), ' MHz'], 'FontSize', fsize, 'Interpreter', 'Latex');
xlabel('Frequency Offset (kHz)', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Phase (degrees)', 'FontSize',fsize, 'Interpreter', 'Latex');
grid on;

d  = 1e3;         % Distance to receiver
phi_freq_error = wrapToPi(2*pi*d.*(fc + df)/c)*180/pi;
plot(1e-3*df, phi_freq_error, 'LineWidth', 2, 'DisplayName', ['d = ', num2str(d), ' m']);

d  = 1e4;         % Distance to receiver
phi_freq_error = wrapToPi(2*pi*d.*(fc + df)/c)*180/pi;
plot(1e-3*df, phi_freq_error, 'LineWidth', 2, 'DisplayName', ['d = ', num2str(d), ' m']);
set(gcf, 'color', 'w');

legend('show', 'FontSize', fsize);

fn = '/Users/Ivan/Documents/Thesis/figures/propagation_phase_vs_freq_5GHz';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');

% -------------------------------------------------------------------------
%------------------------- Position Errors --------------------------------
fc = 915e6;
dx = 0:10:500;
phi_pos_error = wrapToPi(2*pi.*dx*fc/c)*180/pi;

figure();
plot(dx, phi_pos_error, 'LineWidth', 2, 'DisplayName', ['f_c = ', num2str(fc/1e6), ' MHz']);
hold on;
title('Phase vs Position Error', 'FontSize', fsize, 'Interpreter', 'Latex');
xlabel('Position Error (m)', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Phase (degrees)', 'FontSize', fsize, 'Interpreter', 'Latex');
grid on;

fc = 5e9;
phi_pos_error = wrapToPi(2*pi.*dx*fc/c)*180/pi;

plot(dx, phi_pos_error, 'LineWidth', 2, 'DisplayName', ['f_c = ', num2str(fc/1e6), ' MHz']);
legend('show', 'FontSize', fsize);
set(gcf, 'color', 'w');

fn = '/Users/Ivan/Documents/Thesis/figures/propagation_phase_vs_pos';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');