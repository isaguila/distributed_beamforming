classdef ofdm_demodulator
    properties
        fft_len;
        guard_len;
        ofdm_sym_len;
        constellation_syms;
        constellation_bits;
        guard_int;
        ifft_samples;
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
            obj.ifft_samples = obj.syms(obj.guard_len+1:end, :);
            obj.constellation_syms = fftshift(fft(obj.ifft_samples, obj.fft_len))/sqrt(obj.fft_len);
            obj.constellation = obj.constellation.demod_sig(obj.constellation_syms, scale_factor);
            obj.constellation_bits = obj.constellation.bits;
            
            % Parallel to serial conversion
            obj.constellation.bits = reshape(obj.constellation_bits, [], 1);
            obj.constellation.syms = reshape(obj.constellation_syms, [], 1);
        end
    end
end
