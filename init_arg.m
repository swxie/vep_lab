%��ʼ�������ṹ��ʹ�õĺ���
%�������������signal_enhance�����и���
function [arg] = init_arg()
arg = struct('fs', 100, 'delta_t', 1, 'interest_length', 0.3, ...%������Ϣ���ǿɵ�
        'dwt_order', 1, 'wave', "db1", ...%С���任���
        'slice', 1, 'max_iteration', 6, 'threshold', 0.2, 'log_threshold', 3);%�ź���ǿ�����
end

