close all;
clear;
clc;

fsize = 14;

phi_0 = 0;
df = 50:1000;

phi_60 = pi/3;
T_max_60 = 1e3*(phi_60 - phi_0)./(2*pi*df);

phi_30 = pi/6;
T_max_30 = 1e3*(phi_30 - phi_0)./(2*pi*df);

figure();
plot(df, T_max_30, 'LineWidth', 2, 'DisplayName', '$\phi(t) = \pi/6$');
hold on;
plot(df, T_max_60, 'LineWidth', 2, 'DisplayName', '$\phi(t) = \pi/3$');
set(gcf,'color','w');
title('Maximum Packet Length vs CFO', 'FontSize', fsize, 'Interpreter', 'Latex');
xlabel('Frequency Offset (Hz)', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Packet Length (ms)', 'FontSize', fsize, 'Interpreter', 'Latex');
legend('show', 'FontSize', fsize, 'Interpreter', 'Latex');
grid on;

fn = '/Users/Ivan/Documents/Thesis/figures/packet_len_vs_cfo';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');


T_packet = 0.1e-3:1e-8:1.5e-3;
df_max_30 = (phi_30 - phi_0)./(2*pi*T_packet);
df_max_60 = (phi_60 - phi_0)./(2*pi*T_packet);

figure();
plot(1e3*T_packet, df_max_30, 'LineWidth', 2, 'DisplayName', '$\phi(t) = \pi/6$');
hold on;
plot(1e3*T_packet, df_max_60, 'LineWidth', 2, 'DisplayName', '$\phi(t) = \pi/3$');
set(gcf,'color','w');
title('Maximum CFO vs Packet Length', 'FontSize', fsize, 'Interpreter', 'Latex');
xlabel('Packet Length (ms)', 'FontSize', fsize, 'Interpreter', 'Latex');
ylabel('Frequency Offset (Hz)', 'FontSize', fsize, 'Interpreter', 'Latex');
legend('show', 'FontSize', fsize, 'Interpreter', 'Latex');
grid on;

fn = '/Users/Ivan/Documents/Thesis/figures/cfo_vs_packet_len';
[imageData, alpha] = export_fig(fn, '-a1', '-pdf');