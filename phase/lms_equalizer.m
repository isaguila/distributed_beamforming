% function z = lms_equalizer(y, x, num_points, num_taps, mu)
% Avetis Ioannisyan (2020). 
% LMS Adaptation Training Equalizer 
% (https://www.mathworks.com/matlabcentral/fileexchange/9205-lms-adaptation-training-equalizer), 
% MATLAB Central File Exchange. Retrieved May 25, 2020.

function [z, e, w] = lms_equalizer(y, x, num_points, num_taps, mu)

if nargin == 2 % Default parameters
    num_points = length(x);  % # of points
    num_taps = 9;            % channel order
    mu = 0.01;              % iteration step size
end

z = complex(zeros(num_points, 1));
e = complex(zeros(num_points, 1));
w = zeros(num_taps+1,1) + 1i*zeros(num_taps+1,1);

% LMS Adaptation
for n  = num_taps+1 : num_points
    % select part of training input
    in = x(n : -1 : n-num_taps);
    z(n) = w'*in;
    % compute error
    e(n) = y(n)-z(n);
    % update taps
    w = w + mu*(real(e(n)*conj(in)) - 1i*imag(e(n)*conj(in)));
end




end