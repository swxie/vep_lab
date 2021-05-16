%初始化参数结构体使用的函数
%参数的意义会在signal_enhance函数中给出
function [arg] = init_arg()
arg = struct('fs', 100, 'delta_t', 1, 'interest_length', 0.3, ...%传递信息，非可调
        'dwt_order', 1, 'wave', "db1", ...%小波变换相关
        'slice', 1, 'max_iteration', 6, 'threshold', 0.2, 'log_threshold', 3);%信号增强器相关
end

