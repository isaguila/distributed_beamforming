classdef qam_demodulator
    properties
        bits;
        syms;
        M;
    end
    methods
        function obj = qam_demodulator(M)
            obj.M = M;
        end
        function obj = demod_sig(obj, x, scale_factor)
            obj.syms = x*diag(scale_factor./vecnorm(x));
            obj.bits = qamdemod(obj.syms, obj.M, 'OutputType', 'bit', 'UnitAveragePower', true);
        end
    end
end
