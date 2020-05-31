classdef ofdm_modulator
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
        function obj = ofdm_modulator(fft_len, guard_len, constellation)
            obj.fft_len = fft_len;
            obj.guard_len = guard_len;
            obj.ofdm_sym_len = obj.fft_len + obj.guard_len;
            obj.constellation = constellation;
        end
        function obj = get_syms(obj, num_bits)
            obj.constellation = obj.constellation.get_syms(num_bits);
            obj.bits = obj.constellation.bits;
            obj.constellation_syms = obj.constellation.syms;
            obj.constellation_matrix = reshape(obj.constellation_syms, obj.fft_len, length(obj.constellation_syms)/obj.fft_len);
            obj.constellation_matrix = ifft(ifftshift(obj.constellation_matrix), obj.fft_len);
            obj.guard_int = obj.constellation_matrix(end-obj.guard_len+1:end, :);
            obj.syms = [obj.guard_int; obj.constellation_matrix];
        end
    end
end
