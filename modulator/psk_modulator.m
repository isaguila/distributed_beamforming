classdef psk_modulator
    properties
        bits;
        syms;
        M;
    end
    methods
        function obj = psk_modulator(M)
            if (M ~= 2)
                error([num2str(M) '-PSK not implemented yet']);
            end
            obj.M = M;
        end
        function obj = get_syms(obj, num_bits)
            obj.bits = randi([0 1], num_bits, 1);
            switch obj.M
                case 2
                    obj.syms = 2*obj.bits-1;
            end
        end
    end
end
