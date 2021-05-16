function [xcorr_after, xcorr_ori] = testInnerPatient(name, isPlot, arg)
% name示例：'.\data\xsw_?.mat'
% 会根据name读取满足条件的所有mat文件（通常是关于同一个受试者的全部mat文件）
% 每个文件内包含了一个60s×4的四路fvep信号数据，采样率为1e2
if nargin == 1
    isPlot = 1;
    arg = init_arg();
end
filelist = dir(name);
file_num  = size(filelist,1);
left_vep_before = [];
right_vep_before = [];
left_vep_after = [];
right_vep_after = [];
for i = 1:file_num
    %读取文件
    load(strcat('./data/', filelist(i).name), 'data');
    [left_vep_cur, right_vep_cur, learning_curve, left_vep_cur_b, right_vep_cur_b] = signal_enhance(data, arg);
    left_vep_before = [left_vep_before, mean(left_vep_cur_b, 2)];
    right_vep_before = [right_vep_before, mean(right_vep_cur_b, 2)];
    left_vep_after = [left_vep_after, mean(left_vep_cur, 2)];
    right_vep_after = [right_vep_after, mean(right_vep_cur, 2)];
end
if isPlot == 1
    figure;
    subplot(2,1,1);
    plot(mean(left_vep_before,2),'-r');
    hold on;
    plot(mean(left_vep_after,2),'-b');
    legend("before", "after");
    title("left");
    subplot(2,1,2);
    plot(mean(right_vep_before,2),'-r');
    hold on;
    plot(mean(right_vep_after,2),'-b');
    legend("before", "after");
    title("right");
end
xcorr_ori = (mean(mean(corrcoef(left_vep_before))) + mean(mean(corrcoef(right_vep_before)))) / 2;
xcorr_after = (mean(mean(corrcoef(left_vep_after))) + mean(mean(corrcoef(right_vep_after)))) / 2;
end

