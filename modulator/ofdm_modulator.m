classdef ofdm_modulator
    properties
        fft_len;
        guard_len;
        ofdm_sym_len;
        guard_int;
        ifft_samples;
        syms;
        constellation;
        constellation_bits;
        constellation_syms;
    end
    methods
        function obj = ofdm_modulator(fft_len, guard_len, constellation)
            obj.fft_len = fft_len;
            obj.guard_len = guard_len;
            obj.ofdm_sym_len = obj.fft_len + obj.guard_len;
            obj.constellation = constellation;
        end
        function obj = get_syms(obj, num_bits)
            obj.constellation = obj.constellation.get_syms(num_bits);
            num_bits = length(obj.constellation.bits);
            num_constellation_syms = num_bits/log2(obj.constellation.M);
            num_ofdm_syms = num_constellation_syms/obj.fft_len;
            obj.constellation_bits = reshape(obj.constellation.bits, num_bits/num_ofdm_syms, num_ofdm_syms);
            obj.constellation_syms = reshape(obj.constellation.syms, num_constellation_syms/num_ofdm_syms, num_ofdm_syms);
            obj.ifft_samples = ifft(ifftshift(obj.constellation_syms), obj.fft_len)*sqrt(obj.fft_len);
            obj.guard_int = obj.ifft_samples(end-obj.guard_len+1:end, :);
            obj.syms = [obj.guard_int; obj.ifft_samples];
        end
    end
end
