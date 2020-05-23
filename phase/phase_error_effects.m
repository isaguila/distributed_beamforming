close all;
clear;
clc;


%--------------------------------------------------------------------------
%------- Code based on paper by R. Mudumbai, G. Barriac, U. Madhow --------
% On the feasibility of Distributed Beamforming in Wireless Networks
%--------------------------------------------------------------------------
path2save = '/Users/ivan/Documents/Thesis/figures';
epsilon = -180:180;

y = 2*cosd(epsilon/2);

plot(epsilon, y);
title('Combined Signal Amplitude vs Phase Error', 'Interpreter', 'Latex');
xlabel('Phase Error (degrees)', 'Interpreter', 'Latex');
ylabel('Combined Signal Amplitude', 'Interpreter', 'Latex');
grid on;