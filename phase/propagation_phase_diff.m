close all;
clear;
clc;

%--------------------------------------------------------------------------
%------- Code based on paper by R. Mudumbai, G. Barriac, U. Madhow --------
% On the Feasibility of Distributed Beamforming in Wireless Networks
%--------------------------------------------------------------------------
path2save = '/Users/ivan/Documents/Thesis/figures';
Ntx       = 2;          % Number of transmitter nodes
fc        = 915e6;      % carrier frequency in Hz
c         = 299792458;  % speed of light in m/s
d         = 1000;       % Distance to receiver
phi_ideal = wrapToPi(2*pi*d*fc/c)*180/pi; % ideal phase in degrees

%------------------------ Frequency Errors --------------------------------

df = -1e4:500:1e4;
phi_freq_error = wrapToPi(2*pi*d.*(fc + df)/c)*180/pi;

figure();
plot(1e-3*df, phi_freq_error);
title('Phase vs Frequency Offset ', 'Interpreter', 'Latex');
xlabel('Frequency Offset (kHz)', 'Interpreter', 'Latex');
ylabel('Phase (degrees)', 'Interpreter', 'Latex');
grid on;

%------------------------- Position Errors --------------------------------
dx = 0:10:1000;
phi_pos_error = wrapToPi(2*pi.*dx*fc/c)*180/pi;

figure();
plot(dx, phi_pos_error);
title('Phase vs Position Error', 'Interpreter', 'Latex');
xlabel('Position Error (m)', 'Interpreter', 'Latex');
ylabel('Phase (degrees)', 'Interpreter', 'Latex');
grid on;