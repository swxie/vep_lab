%��ʼ�������ṹ��ʹ�õĺ���
%�������������signal_enhance�����и���
function [arg] = init_arg()
arg = struct('fs', 100, 'delta_t', 1, 'interest_length', 0.3, ...%������Ϣ���ǿɵ�
        'dwt_order', 1, 'wave', "db1", ...
        'slice', 2);
end

