classdef psk_demodulator
    properties
        bits;
        syms;
        M;
    end
    methods
        function obj = psk_demodulator(M)
            if (M ~= 2)
                error([num2str(M) '-PSK not implemented yet']);
            end
            obj.M = M;
        end
        function obj = demod_sig(obj, x, scale_factor)
            obj.syms = x*diag(scale_factor/vecnorm(x));
            zeros_idx = real(obj.syms) <= 0;
            ones_idx = real(obj.syms) > 0;
            obj.bits = zeros(size(obj.syms));
            obj.bits(zeros_idx) = 0;
            obj.bits(ones_idx) = 1;
        end
    end
end
