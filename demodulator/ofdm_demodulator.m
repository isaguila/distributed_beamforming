classdef ofdm_demodulator
    properties
        fft_len;
        guard_len;
        ofdm_sym_len;
        bits;
        constellation_syms;
        constellation_matrix;
        guard_int;
        syms;
        constellation;
    end
    methods
        function obj = ofdm_demodulator(fft_len, guard_len, constellation)
            obj.fft_len = fft_len;
            obj.guard_len = guard_len;
            obj.ofdm_sym_len = obj.fft_len + obj.guard_len;
            obj.constellation = constellation;
        end
        function obj = demod_sig(obj, x, scale_factor)
            sig_len = length(x);
            obj.syms = reshape(x, obj.ofdm_sym_len, sig_len/obj.ofdm_sym_len);
            obj.constellation_matrix = obj.syms(obj.guard_len+1:end, :);
            obj.constellation_matrix = fftshift(fft(obj.constellation_matrix, obj.fft_len));
            obj.constellation_syms = reshape(obj.constellation_matrix, [], 1);
            obj.constellation_syms = obj.constellation_syms*scale_factor/norm(obj.constellation_syms);
            obj.bits = qamdemod(obj.constellation_syms, obj.constellation.M, 'OutputType', 'bit', 'UnitAveragePower', true);
        end
    end
end
