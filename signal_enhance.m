function [left_vep, right_vep, learning_curve, left_vep_before, right_vep_before] = signal_enhance(data, arg)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
%输入vep四路信号输入，返回合成的left_vep和right_vep
%   定义全局变量
if nargin == 1
    arg = init_arg();
end
%设置参数
fs  = arg.fs;%频率
delta_t = arg.delta_t;%刺激的时间间隔
window_length = delta_t * fs;%刺激间隔对应的窗口长度
interest_length = arg.interest_length * fs;%感兴趣区间长度
%设置小波变换参数
dwt_order = arg.dwt_order;%dwt阶数
wave = arg.wave;%小波类型
%筛选相关参数
slice = arg.slice;%每次筛选筛选出的低质量信号片段个数
max_iteration = arg.max_iteration;%最大迭代次数
threshold = arg.threshold;%迭代终止对应的下降值
log_threshold = arg.log_threshold;%迭代终止时的信号功率比
%构造分段
buffer = data(1:window_length, 1);
diff = conv(buffer,[1;-1]);
diff(1:5) = zeros(1,5);
diff(end - 4 : end) = zeros(1,5);
max_index = find(diff == max(diff),1);
end_mat = (max_index + interest_length - 1) : window_length : size(data,1);
begin_mat = end_mat - interest_length + 1;
left_vep = zeros(interest_length, length(begin_mat));
right_vep = zeros(interest_length, length(begin_mat));
%构造分类样本
for i = 1:length(begin_mat)
    left_vep(:, i) = data(begin_mat(i) : end_mat(i), 2);
    right_vep(:, i) = data(begin_mat(i) : end_mat(i), 3);
end
left_vep_before = left_vep;
right_vep_before = right_vep;
%%计算dwt后长度
temp = left_vep(:,1);
for k = 1 : dwt_order
    temp = dwt(temp,wave);
end
dwt_length = length(temp);
%%初始化vep矩阵
vep = zeros(2 * dwt_length, length(begin_mat));
for i = 1:length(begin_mat)
    left_temp = left_vep(:, i);
    right_temp = right_vep(:, i);
    for k = 1 : dwt_order
        left_temp = dwt(left_temp,wave);
        right_temp = dwt(right_temp, wave);
    end
    vep(:, i) = [left_temp; right_temp];
end
% figure;
% subplot(2, 2, 1);
% plot(sum(left_vep, 2));
% subplot(2, 2, 2);
% plot(sum(right_vep, 2));
%%进行迭代判别、删除
learning_curve = [];
for i = 1 : max_iteration
    mask = svm_classify(vep, slice);
    left_vep_new = mean(left_vep(:, mask > 0), 2);
    right_vep_new = mean(right_vep(:, mask > 0), 2);
    left_vep_del = mean(left_vep(:, mask <= 0), 2);
    right_vep_del = mean(right_vep(:, mask <= 0), 2);
    %条件一求去除后的互相关
    left_vep_target = [corrcoef(left_vep_new, left_vep_del)];
    right_vep_target = [corrcoef(right_vep_new, right_vep_del)];
    left_vep_target = left_vep_target(2, 1);
    right_vep_target = right_vep_target(2, 1);
    %条件2：两个信号功率之差不应该过大
    left_ratio = log(mean(left_vep_new.^2) / mean(left_vep_del .^2));
    right_ratio = log(mean(right_vep_new.^2) / mean(right_vep_del .^2));
    learning_curve = [learning_curve, min([left_vep_target, right_vep_target])];
    if learning_curve(end) > threshold && max(abs(left_ratio), abs(right_ratio)) < log_threshold
        break;
    end
    left_vep = left_vep(:, mask > 0);
    right_vep = right_vep(:, mask > 0);
    vep = vep(:, mask > 0);
end
end

function [mask] = svm_classify(data, slice)
    SVMModel = fitcsvm(data', ones(size(data, 2), 1),'kernelScale',0.04, ...
            'OutlierFraction',0.05, 'nu', 0.5);
    [~,scorePred] = predict(SVMModel, data');
    t = sort(scorePred);
    if length(scorePred) < slice + 1
        mask = ones(size(scorePred));
    else
        mask = scorePred >= t(slice + 1);
    end
end

