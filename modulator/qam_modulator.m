classdef qam_modulator
    properties
        bits;
        syms;
        M;
    end
    methods
        function obj = qam_modulator(M)
            obj.M = M;
        end
        function obj = get_syms(obj, num_bits)
            obj.bits = randi([0 1], num_bits, 1);
            switch obj.M
                case 4
                    obj.syms = qammod(obj.bits, obj.M, 'InputType', 'bit', 'UnitAveragePower', true);
                case 16
                    obj.syms = qammod(obj.bits, obj.M, 'InputType', 'bit', 'UnitAveragePower', true);
                case 64
                    obj.syms = qammod(obj.bits, obj.M, 'InputType', 'bit', 'UnitAveragePower', true);
                case 256
                    obj.syms = qammod(obj.bits, obj.M, 'InputType', 'bit', 'UnitAveragePower', true);
                case 1024
                    obj.syms = qammod(obj.bits, obj.M, 'InputType', 'bit', 'UnitAveragePower', true);
            end
        end
    end
end
