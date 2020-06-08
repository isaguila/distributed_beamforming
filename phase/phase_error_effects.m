close all;
clear;
clc;

%--------------------------------------------------------------------------
%------- Code based on paper by R. Mudumbai, G. Barriac, U. Madhow --------
% On the feasibility of Distributed Beamforming in Wireless Networks
%--------------------------------------------------------------------------
fn = '/Users/ivan/Documents/Thesis/figures/phase_error_effects';
epsilon = -180:180;

fsize = 14;

y = 2*cosd(epsilon/2);

plot(epsilon, y, 'LineWidth', 2);
title('Combined Signal Amplitude vs Phase Error', 'FontSize', fsize, 'Interpreter', 'Latex');
xlabel('Phase Error (degrees)', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Combined Signal Amplitude', 'FontSize', fsize, 'Interpreter', 'Latex');
grid on;
set(gcf, 'color', 'w');

[imageData, alpha] = export_fig(fn, '-a1', '-pdf');


