function evm = get_evm_qam(x, M, scale)
% Non genie-aided EVM calculation based on nearest decision region

if (nargin < 3)
    scale = 'dB';
end

% Generate ideal M-ary QAM symbols
ideal_syms = qammod((0:M-1).', M, 'UnitAveragePower', true);

% Repeat the received symbols
x_matrix = repmat(x.', length(ideal_syms), 1);

% Repeat ideal symbols to match dimensions of received signal
ideal_matrix = repmat(ideal_syms, 1, length(x));

% Calculate minimum Euclidean distance between received symbol and closest
% decision region
sym_diff = min(abs(x_matrix - ideal_matrix));

% Calculate the EVM
switch scale
    case 'dB'
        evm = 10*log10(mean(sym_diff.^2));
    case 'Linear'
        evm = mean(sym_diff.^2);
end